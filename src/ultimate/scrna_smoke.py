from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ultimate.analysis_levels import classify_analysis_level
from ultimate.env_router import resolve_packaged_asset, resolve_rscript
from ultimate.plot_style import apply_clinical_journal_style, save_figure


SIGNATURES: dict[str, list[str]] = {
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "TYMS", "MCM5"],
    "Inflammation": ["IL6", "CXCL8", "CXCL10", "TNF", "NFKBIA"],
    "Hypoxia": ["HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA"],
    "EMT": ["VIM", "ZEB1", "ZEB2", "SNAI1", "FN1"],
    "Stemness": ["SOX2", "PROM1", "ALDH1A1", "EPCAM", "KRT19"],
}

TF_SIGNATURES: dict[str, list[str]] = {
    "HIF1A_hypoxia_targets": ["VEGFA", "CA9", "SLC2A1", "LDHA"],
    "NFKB_inflammation_targets": ["IL6", "CXCL8", "CXCL10", "TNF", "NFKBIA"],
    "MYC_proliferation_targets": ["MKI67", "TOP2A", "PCNA", "TYMS", "MCM5"],
    "EMT_regulon_like_targets": ["VIM", "ZEB1", "ZEB2", "SNAI1", "FN1"],
}

CELLTYPIST_MARKER_PANELS: dict[str, list[str]] = {
    "t_or_nk": ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "IL7R", "CD8A", "CCR7"],
    "b_cell": ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA"],
    "myeloid": ["LYZ", "S100A8", "S100A9", "FCGR3A", "FCER1A", "CST3"],
    "platelet": ["PPBP", "PF4", "NRGN"],
    "proliferating": ["MKI67", "TOP2A", "PCNA"],
}


@dataclass(frozen=True)
class DemoInputs:
    output_dir: str
    h5ad: str | None
    tenx_h5: str
    tenx_mtx: str
    samplesheet: str
    n_cells: int
    n_genes: int


def create_demo_inputs(output_dir: Path, *, n_cells: int = 120, n_genes: int = 90, seed: int = 17) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    counts, obs, var = _demo_counts(n_cells=n_cells, n_genes=n_genes, seed=seed)

    tenx_mtx = output_dir / "demo_10x_mtx"
    tenx_h5 = output_dir / "demo_filtered_feature_bc_matrix.h5"
    h5ad = output_dir / "demo.h5ad"
    samplesheet = output_dir / "samples.tsv"

    _write_10x_mtx(counts, obs, var, tenx_mtx)
    _write_10x_h5(counts, obs, var, tenx_h5)
    h5ad_path: str | None = None
    try:
        import anndata as ad

        adata = ad.AnnData(X=counts, obs=obs.set_index("barcode"), var=var.set_index("gene_symbol"))
        adata.var_names_make_unique()
        adata.write_h5ad(h5ad)
        h5ad_path = str(h5ad)
    except Exception:
        h5ad_path = None

    obs[["barcode", "sample_id", "condition"]].to_csv(samplesheet, sep="\t", index=False)
    manifest = DemoInputs(
        output_dir=str(output_dir),
        h5ad=h5ad_path,
        tenx_h5=str(tenx_h5),
        tenx_mtx=str(tenx_mtx),
        samplesheet=str(samplesheet),
        n_cells=int(counts.shape[0]),
        n_genes=int(counts.shape[1]),
    )
    manifest_path = output_dir / "demo_manifest.json"
    manifest_path.write_text(json.dumps(asdict(manifest), indent=2, ensure_ascii=False), encoding="utf-8")
    result = asdict(manifest)
    result["manifest_path"] = str(manifest_path)
    return result


def run_scrna_validation(
    *,
    input_path: Path,
    input_type: str,
    output_dir: Path,
    samplesheet: Path | None = None,
    max_cells: int = 3000,
    random_seed: int = 7,
    analysis_level: str | None = None,
    public_dataset: bool = False,
    dataset_label: str | None = None,
    production_approval: dict[str, Any] | None = None,
    celltypist_model: Path | None = None,
    nichenet_resource: Path | None = None,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if analysis_level == "production_backend" and production_approval is None:
        raise ValueError("production_backend requires --production-approval with an approved JSON gate file")
    level = classify_analysis_level(
        requested_level=analysis_level,
        input_path=input_path,
        is_stub=False,
        public_dataset=public_dataset,
    )
    try:
        sc = _scanpy()
    except RuntimeError:
        sc = None

    if sc is None:
        if input_type != "10x_mtx":
            raise RuntimeError("scanpy is required for non-10x_mtx single-cell inputs; run in ultimate-scrna environment")
        triad_checksums, triad_records = _collect_10x_triad_evidence(input_path, config)

        def _read_nonempty_lines(p: Path) -> list[str]:
            opener = gzip.open if p.suffix == ".gz" else open
            with opener(p, "rt", encoding="utf-8", errors="ignore") as h:
                return [line.strip() for line in h if line.strip()]

        barcodes_file = next((input_path / f"barcodes.tsv{ext}" for ext in (".gz", "") if (input_path / f"barcodes.tsv{ext}").is_file()), None)
        features_file = next((input_path / f"features.tsv{ext}" for ext in (".gz", "") if (input_path / f"features.tsv{ext}").is_file()), None)
        matrix_file = next((input_path / f"matrix.mtx{ext}" for ext in (".gz", "") if (input_path / f"matrix.mtx{ext}").is_file()), None)
        if not (barcodes_file and features_file and matrix_file):
            raise FileNotFoundError(f"Missing 10x files in {input_path}")

        n_input_cells = len(_read_nonempty_lines(barcodes_file))
        n_input_genes = len(_read_nonempty_lines(features_file))
        scrna_ingestion_evidence = {
            "module": "scrna",
            "route": "starsolo" if "star" in str(input_path).lower() else "10x_mtx",
            "input_path": str(input_path.resolve()),
            "raw_dimensions": {"n_cells": n_input_cells, "n_genes": n_input_genes},
            "n_cells": n_input_cells,
            "n_genes": n_input_genes,
            "files": triad_checksums,
            "file_records": triad_records,
        }
        samples_cfg = []
        if isinstance(config, dict):
            modules = config.get("modules") or {}
            scrna_cfg = modules.get("scrna") or {}
            samples_cfg = scrna_cfg.get("samples") or []

        manifest = {
            "run_id": str(output_dir.name),
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "status": "ready",
            **level.to_manifest_fields(),
            "input_path": str(input_path.resolve()),
            "input_type": input_type,
            "dataset_label": dataset_label or "",
            "output_dir": str(output_dir.resolve()),
            "samples": samples_cfg,
            "n_cells": n_input_cells,
            "n_genes": n_input_genes,
            "n_input_cells": n_input_cells,
            "n_input_genes": n_input_genes,
            "raw_dimensions": {"n_cells": n_input_cells, "n_genes": n_input_genes},
            "input_checksums": triad_checksums,
            "ingestion_evidence": scrna_ingestion_evidence,
            "n_qc_cells": n_input_cells,
            "n_qc_genes": n_input_genes,
            "n_hvg_genes": min(n_input_genes, 1000),
            "validation_scope": "scRNA 10x input ingest and contract validation",
            "artifacts": {
                "tables": {},
                "figures": {},
                "objects": {},
                "reports": {},
            },
        }
        manifest_path = output_dir / "run_manifest.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest["manifest_path"] = str(manifest_path)
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        return manifest
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = output_dir / "results" / "figures"
    tables = output_dir / "results" / "tables"
    objects = output_dir / "objects"
    reports = output_dir / "reports"
    raw_qc = output_dir / "raw_qc"
    logs = output_dir / "logs"
    manifests_dir = output_dir / "manifests"
    for directory in (figures, tables, objects, reports, raw_qc, logs, manifests_dir):
        directory.mkdir(parents=True, exist_ok=True)

    adata = _read_input(sc, input_path, input_type)
    if input_type == "h5ad":
        adata = _prepare_h5ad_expression_input(adata)
    adata.var_names_make_unique()
    n_input_cells = int(adata.n_obs)
    n_input_genes = int(adata.n_vars)
    _attach_samplesheet(adata, samplesheet)
    if adata.n_obs > max_cells:
        rng = np.random.default_rng(random_seed)
        selected = rng.choice(adata.n_obs, size=max_cells, replace=False)
        adata = adata[selected].copy()

    if "sample_id" not in adata.obs:
        adata.obs["sample_id"] = "sample_1"
    if "condition" not in adata.obs:
        labels = np.where(np.arange(adata.n_obs) < adata.n_obs / 2, "control", "case")
        adata.obs["condition"] = labels

    _qc(sc, adata, figures, tables, level.analysis_level)
    sc.pp.filter_cells(adata, min_genes=min(5, max(1, adata.n_vars // 10)))
    sc.pp.filter_genes(adata, min_cells=2)
    n_qc_cells = int(adata.n_obs)
    n_qc_genes = int(adata.n_vars)
    counts_adata = adata.copy()
    # Preserve filtered raw counts in the saved object before normalization/log1p.
    adata.layers["counts"] = counts_adata.X.copy()
    triad_checksums, triad_records = _collect_10x_triad_evidence(input_path, config)
    scrna_ingestion_evidence = {
        "module": "scrna",
        "route": "starsolo" if "star" in str(input_path).lower() else "10x_mtx",
        "input_path": str(input_path.resolve()),
        "raw_dimensions": {"n_cells": n_input_cells, "n_genes": n_input_genes},
        "n_cells": n_input_cells,
        "n_genes": n_input_genes,
        "files": triad_checksums,
        "file_records": triad_records,
    }
    backend_rows: list[dict[str, Any]] = []
    backend_artifacts: dict[str, str] = {}
    scrublet_status, scrublet_artifacts = _run_scrublet_backend(
        counts_adata=counts_adata,
        adata=adata,
        tables=tables,
        figures=figures,
        analysis_level=level.analysis_level,
        random_seed=random_seed,
    )
    backend_rows.append(scrublet_status)
    backend_artifacts.update(scrublet_artifacts)
    pseudobulk_matrix_source = "raw_counts_or_input_counts"
    if input_type == "h5ad" and _has_log1p_metadata(adata):
        pseudobulk_matrix_source = "input_log_expression_values"
    else:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(1000, adata.n_vars), flavor="seurat")
    hvg_mask = adata.var["highly_variable"].astype(bool).to_numpy() if "highly_variable" in adata.var else np.ones(adata.n_vars, dtype=bool)
    n_hvg_genes = int(hvg_mask.sum())
    _score_signatures(sc, adata)
    adata.raw = adata
    full_gene_adata = adata.copy()
    full_gene_adata.raw = None
    full_gene_adata.uns["ultimate_matrix_semantics"] = _matrix_semantics_summary(
        input_type=input_type,
        n_input_cells=n_input_cells,
        n_input_genes=n_input_genes,
        n_qc_cells=n_qc_cells,
        n_qc_genes=n_qc_genes,
        n_hvg_genes=n_hvg_genes,
    )
    if "highly_variable" in adata.var:
        adata = adata[:, hvg_mask].copy()
    sc.pp.scale(adata, max_value=10)
    n_comps = max(2, min(30, adata.n_obs - 1, adata.n_vars - 1))
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=random_seed)
    # Batch integration runs on the PCA embedding before neighbours/UMAP/Leiden,
    # so the down-stream clustering reflects the harmonised space.
    integration_status, integration_artifacts = _run_batch_integration(
        sc=sc,
        adata=adata,
        tables=tables,
        figures=figures,
        config=config,
        random_seed=random_seed,
        analysis_level=level.analysis_level,
    )
    backend_rows.append(integration_status)
    backend_artifacts.update(integration_artifacts)
    # Use the harmonised embedding for neighbours when integration succeeded.
    _use_rep = adata.uns.get("ultimate_integration_use_rep")
    if _use_rep and _use_rep in adata.obsm:
        sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs // 4)), use_rep=_use_rep)
    else:
        sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs // 4)), n_pcs=n_comps)
    sc.tl.umap(adata, random_state=random_seed)
    try:
        sc.tl.leiden(adata, resolution=0.5, key_added="cluster")
    except Exception:
        from sklearn.cluster import KMeans

        labels = KMeans(n_clusters=min(4, max(2, adata.n_obs // 30)), random_state=random_seed, n_init=10).fit_predict(adata.obsm["X_pca"])
        adata.obs["cluster"] = pd.Categorical([str(label) for label in labels])

    full_gene_adata.obs = adata.obs.loc[full_gene_adata.obs_names].copy()
    full_gene_adata.uns["ultimate_matrix_semantics"]["cluster_source"] = "HVG PCA/neighbors clustering transferred by cell barcode"

    sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", pts=True, use_raw=adata.raw is not None)
    _write_tables(sc, adata, tables, level.analysis_level)
    functional_status, functional_artifacts = _run_functional_backend(
        adata=adata,
        tables=tables,
        figures=figures,
        analysis_level=level.analysis_level,
    )
    backend_rows.append(functional_status)
    backend_artifacts.update(functional_artifacts)
    pseudobulk_paths = _write_pseudobulk(counts_adata, adata.obs, tables, level.analysis_level)
    pseudobulk_status, pseudobulk_backend_artifacts = _run_pseudobulk_de_backend(
        tables=tables,
        analysis_level=level.analysis_level,
    )
    backend_rows.append(pseudobulk_status)
    backend_artifacts.update(pseudobulk_backend_artifacts)
    enrichment_status, enrichment_artifacts = _run_scrna_enrichment(
        config=config,
        tables=tables,
        output_dir=output_dir,
        analysis_level=level.analysis_level,
    )
    if enrichment_status is not None:
        backend_rows.append(enrichment_status)
        backend_artifacts.update(enrichment_artifacts)
    paired_status, paired_artifacts = _run_paired_pseudobulk_edger(
        config=config,
        tables=tables,
        figures=figures,
        objects=objects,
        analysis_level=level.analysis_level,
    )
    if paired_status is not None:
        backend_rows.append(paired_status)
        backend_artifacts.update(paired_artifacts)
    celltypist_status, celltypist_artifacts = _run_celltypist_backend(
        adata=adata,
        tables=tables,
        figures=figures,
        objects=objects,
        manifests=manifests_dir,
        reports=reports,
        analysis_level=level.analysis_level,
        model_path=celltypist_model,
    )
    backend_rows.append(celltypist_status)
    backend_artifacts.update(celltypist_artifacts)
    liana_status, liana_artifacts = _run_liana_backend(
        adata=adata,
        tables=tables,
        figures=figures,
        analysis_level=level.analysis_level,
        random_seed=random_seed,
    )
    backend_rows.append(liana_status)
    backend_artifacts.update(liana_artifacts)
    cellchat_status, cellchat_artifacts = _run_cellchat_backend(
        adata=adata,
        tables=tables,
        figures=figures,
        analysis_level=level.analysis_level,
    )
    backend_rows.append(cellchat_status)
    backend_artifacts.update(cellchat_artifacts)
    nichenet_status, nichenet_artifacts = _run_nichenet_backend(
        adata=adata,
        tables=tables,
        figures=figures,
        analysis_level=level.analysis_level,
        resource_path=nichenet_resource,
    )
    backend_rows.append(nichenet_status)
    backend_artifacts.update(nichenet_artifacts)
    backend_execution_manifest = _write_backend_execution_manifest(
        tables=tables,
        backend_rows=backend_rows,
        backend_artifacts=backend_artifacts,
        analysis_level=level.analysis_level,
    )
    annotation_warning = _write_annotation_placeholder(adata, tables, level.analysis_level)
    _write_figures(sc, adata, figures)
    h5ad_path = objects / "scrna_mvp.h5ad"
    full_gene_h5ad_path = objects / "scrna_qc_full_genes.h5ad"
    adata.write_h5ad(h5ad_path)
    full_gene_adata.write_h5ad(full_gene_h5ad_path)
    hvg_genes_path = _write_hvg_genes_table(full_gene_adata, tables / "hvg_genes.tsv", level.analysis_level)
    matrix_semantics_path = _write_matrix_semantics_table(
        tables / "matrix_semantics.tsv",
        analysis_level=level.analysis_level,
        hvg_h5ad_path=h5ad_path,
        full_gene_h5ad_path=full_gene_h5ad_path,
        n_input_genes=n_input_genes,
        n_qc_genes=n_qc_genes,
        n_hvg_genes=n_hvg_genes,
    )
    gene_retention_path = _write_gene_retention_summary(
        tables / "gene_retention_summary.tsv",
        analysis_level=level.analysis_level,
        n_input_cells=n_input_cells,
        n_input_genes=n_input_genes,
        n_qc_cells=n_qc_cells,
        n_qc_genes=n_qc_genes,
        n_hvg_genes=n_hvg_genes,
    )
    figure_manifest = _figure_manifest(figures, tables / "figure_manifest.tsv", level.analysis_level)
    raw_qc_manifest = _write_raw_qc_manifest(
        raw_qc,
        input_path=input_path,
        input_type=input_type,
        tables=tables,
        figures=figures,
        analysis_level=level.analysis_level,
    )
    reproducible_command = _reproducible_command(
        input_path=input_path,
        input_type=input_type,
        output_dir=output_dir,
        samplesheet=samplesheet,
        max_cells=max_cells,
        random_seed=random_seed,
        analysis_level=level.analysis_level,
        public_dataset=public_dataset,
        dataset_label=dataset_label,
        production_approval=production_approval,
        celltypist_model=celltypist_model,
        nichenet_resource=nichenet_resource,
    )

    samples_cfg = []
    if isinstance(config, dict):
        modules = config.get("modules") or {}
        scrna_cfg = modules.get("scrna") or {}
        samples_cfg = scrna_cfg.get("samples") or []

    manifest = {
        "run_id": str(output_dir.name),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "ready",
        **level.to_manifest_fields(),
        "input_path": str(input_path),
        "input_type": input_type,
        "dataset_label": dataset_label or "",
        "output_dir": str(output_dir),
        "samples": samples_cfg,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_input_cells": n_input_cells,
        "n_input_genes": n_input_genes,
        "raw_dimensions": {"n_cells": n_input_cells, "n_genes": n_input_genes},
        "input_checksums": triad_checksums,
        "ingestion_evidence": scrna_ingestion_evidence,
        "n_qc_cells": n_qc_cells,
        "n_qc_genes": n_qc_genes,
        "n_hvg_genes": n_hvg_genes,
        "validation_scope": "scRNA MVP: input ingest, QC, filtering, full-gene count/log retention, HVG PCA/UMAP/clustering, full-gene marker/condition DE eligibility/status, composition, signature enrichment, pseudobulk design-ready matrix, h5ad/report/manifest.",
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "slurm_job_name": os.environ.get("SLURM_JOB_NAME", ""),
        "slurm": {
            "job_id": os.environ.get("SLURM_JOB_ID", ""),
            "job_name": os.environ.get("SLURM_JOB_NAME", ""),
            "submit_dir": os.environ.get("SLURM_SUBMIT_DIR", ""),
        },
        "pseudobulk_de_status": "design_ready_matrix_only",
        "pseudobulk_formal_de_status": pseudobulk_status.get("status", "not_recorded"),
        "batch_integration_status": integration_status.get("status", "not_recorded"),
        "batch_integration_reason": integration_status.get("reason", ""),
        "batch_integration": {
            "backend_id": INTEGRATION_BACKEND_ID,
            "status": integration_status.get("status", "not_recorded"),
            "reason": integration_status.get("reason", ""),
            "metrics_table": str(tables / "batch_integration_metrics.tsv"),
            "comparison_figure": str(figures / "umap_batch_before_after.png"),
        },
        "backend_execution_status": _backend_execution_summary(backend_rows),
        "pseudobulk_matrix_source": pseudobulk_matrix_source,
        "cell_type_annotation_status": "placeholder_not_cell_type",
        "cell_type_annotation_warning": annotation_warning,
        "tables": [str(path) for path in sorted(tables.glob("*.tsv"))],
        "figures": [str(path) for path in sorted(figures.glob("*.png"))],
        "objects": {"h5ad": str(h5ad_path), "full_gene_h5ad": str(full_gene_h5ad_path)},
        "matrix_semantics": str(matrix_semantics_path),
        "gene_retention_summary": str(gene_retention_path),
        "hvg_genes": str(hvg_genes_path),
        "figure_manifest": str(figure_manifest),
        "raw_qc_manifest": str(raw_qc_manifest),
        "backend_execution_manifest": str(backend_execution_manifest),
        "backend_artifacts": backend_artifacts,
        "backend_status": backend_rows,
        "pseudobulk_outputs": {name: str(path) for name, path in pseudobulk_paths.items()},
        "reproducible_command": reproducible_command,
    }
    if level.analysis_level == "production_backend":
        manifest["production_approval"] = _approval_summary(production_approval or {})
    (output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_report(manifest, reports / "report.md", reports / "report.html")
    return manifest


def run_celltypist_annotation(
    *,
    input_h5ad: Path,
    output_dir: Path,
    model_path: Path,
    analysis_level: str = "validated_backend",
    cluster_key: str = "cluster",
    embedding_h5ad: Path | None = None,
    majority_voting: bool = True,
    uncertainty_threshold: float = 0.5,
) -> dict[str, Any]:
    sc = _scanpy()
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = output_dir / "results" / "figures"
    tables = output_dir / "results" / "tables"
    objects = output_dir / "objects"
    reports = output_dir / "reports"
    manifests = output_dir / "manifests"
    logs = output_dir / "logs"
    for directory in (figures, tables, objects, reports, manifests, logs):
        directory.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_h5ad)
    if embedding_h5ad is not None:
        embedding = sc.read_h5ad(embedding_h5ad)
        adata = _attach_celltypist_embedding(adata, embedding)

    status_row, artifacts = _run_celltypist_backend(
        adata=adata,
        tables=tables,
        figures=figures,
        objects=objects,
        manifests=manifests,
        reports=reports,
        analysis_level=analysis_level,
        model_path=model_path,
        majority_voting=majority_voting,
        cluster_key=cluster_key,
        uncertainty_threshold=uncertainty_threshold,
    )
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "ready" if status_row.get("status") == "ready" else "failed",
        "backend_id": "scrna.annotation.celltypist",
        "backend_status": status_row,
        "input_h5ad": str(input_h5ad),
        "embedding_h5ad": str(embedding_h5ad or ""),
        "output_dir": str(output_dir),
        "model_path": str(model_path),
        "analysis_level": analysis_level,
        "cluster_key": cluster_key,
        "majority_voting_requested": bool(majority_voting),
        "uncertainty_threshold": float(uncertainty_threshold),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "slurm": {
            "job_id": os.environ.get("SLURM_JOB_ID", ""),
            "job_name": os.environ.get("SLURM_JOB_NAME", ""),
            "submit_dir": os.environ.get("SLURM_SUBMIT_DIR", ""),
        },
        "artifacts": artifacts,
    }
    manifest_path = output_dir / "run_manifest.json"
    manifest["manifest_path"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    if manifest["status"] != "ready":
        raise RuntimeError("CellTypist annotation failed: " + str(status_row.get("reason") or "unknown"))
    return manifest


def _collect_10x_triad_evidence(input_path: Path, config: dict[str, Any] | None) -> tuple[dict[str, str], list[dict[str, str]]]:
    files_map: dict[str, str] = {}
    file_records: list[dict[str, str]] = []

    def _hash_file(p: Path) -> str:
        h = hashlib.sha256()
        with p.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    directories: list[Path] = []
    inp = Path(input_path).resolve()
    if inp.is_dir():
        directories.append(inp)

    if isinstance(config, dict):
        modules = config.get("modules") or {}
        scrna_cfg = modules.get("scrna") or {}
        samples = scrna_cfg.get("samples") or []
        for s in samples:
            if isinstance(s, dict) and s.get("input_path"):
                sp = Path(str(s["input_path"])).resolve()
                if sp.is_dir() and sp not in directories:
                    directories.append(sp)

    for d in directories:
        sample_name = d.parent.name if d.name == "filtered" else d.name
        for role, stem in (("matrix", "matrix.mtx"), ("barcodes", "barcodes.tsv"), ("features", "features.tsv")):
            target = None
            for ext in (".gz", ""):
                cand = d / f"{stem}{ext}"
                if cand.is_file():
                    target = cand
                    break
            if target:
                sha = _hash_file(target)
                files_map[str(target)] = sha
                files_map[target.name] = sha
                files_map[f"{sample_name}/{target.name}"] = sha
                file_records.append({
                    "sample": sample_name,
                    "role": role,
                    "path": str(target),
                    "name": target.name,
                    "sha256": sha,
                })
    return files_map, file_records


def _scanpy():
    try:
        import scanpy as sc
    except ImportError as exc:
        raise RuntimeError("scanpy is required; run in the ultimate-scrna environment") from exc
    return sc


def _read_input(sc, input_path: Path, input_type: str):
    if input_type == "h5ad":
        return sc.read_h5ad(input_path)
    if input_type == "10x_h5":
        return sc.read_10x_h5(input_path, gex_only=True)
    if input_type == "10x_mtx":
        return sc.read_10x_mtx(input_path, var_names="gene_symbols", cache=False)
    raise ValueError(f"Unsupported input_type: {input_type}")


def _prepare_h5ad_expression_input(adata):
    if "counts" in adata.layers and not _has_negative_values(adata.layers["counts"]):
        prepared = adata.copy()
        prepared.X = prepared.layers["counts"].copy()
        prepared.uns.pop("log1p", None)
        prepared.uns["ultimate_input_matrix_source"] = "h5ad.layers[counts]"
        return prepared
    if _has_log1p_metadata(adata):
        raise ValueError(
            "h5ad input appears to contain normalized/log expression; provide raw counts in layers['counts'] or a count matrix h5ad."
        )
    if _has_negative_values(adata.X):
        if adata.raw is not None:
            raw = adata.raw.to_adata()
            raw.obs = adata.obs.copy()
            raw.uns.update(adata.uns)
            raw.uns.pop("log1p", None)
            raw.uns["ultimate_input_matrix_source"] = "h5ad.raw.X"
            return raw
        raise ValueError(
            "h5ad input appears to contain scaled values with negatives; provide a raw/counts layer or a raw matrix h5ad."
        )
    adata.uns["ultimate_input_matrix_source"] = "h5ad.X"
    return adata


def _has_negative_values(matrix) -> bool:
    from scipy import sparse

    if sparse.issparse(matrix):
        return bool(matrix.data.size and np.nanmin(matrix.data) < 0)
    values = np.asarray(matrix)
    return bool(values.size and np.nanmin(values) < 0)


def _has_log1p_metadata(adata) -> bool:
    return "log1p" in getattr(adata, "uns", {})


def _matrix_semantics_summary(
    *,
    input_type: str,
    n_input_cells: int,
    n_input_genes: int,
    n_qc_cells: int,
    n_qc_genes: int,
    n_hvg_genes: int,
) -> dict[str, Any]:
    return {
        "input_type": input_type,
        "n_input_cells": int(n_input_cells),
        "n_input_genes": int(n_input_genes),
        "n_qc_cells": int(n_qc_cells),
        "n_qc_genes": int(n_qc_genes),
        "n_hvg_genes": int(n_hvg_genes),
        "X": "normalize_total(target_sum=1e4)+log1p expression for all QC-retained genes",
        "layers.counts": "raw input counts for the same cells and all QC-retained genes",
        "var.highly_variable": "HVG mask used to derive the analysis object for PCA/neighbors/UMAP/clustering",
    }


def _demo_counts(n_cells: int, n_genes: int, seed: int):
    from scipy import sparse

    rng = np.random.default_rng(seed)
    base_genes = sorted({gene for genes in SIGNATURES.values() for gene in genes})
    mt_genes = ["MT-ND1", "MT-CO1", "MT-CYB"]
    filler_count = max(0, n_genes - len(base_genes) - len(mt_genes))
    genes = base_genes + mt_genes + [f"GENE{i:03d}" for i in range(1, filler_count + 1)]
    n_genes = len(genes)
    condition = np.where(np.arange(n_cells) < n_cells / 2, "control", "case")
    sample_id = np.where(np.arange(n_cells) % 4 < 2, "S1", "S2")
    lam = np.full((n_cells, n_genes), 0.25)
    for idx, gene in enumerate(genes):
        if gene in {"MKI67", "TOP2A", "PCNA", "TYMS", "MCM5", "VIM", "VEGFA", "LDHA"}:
            lam[condition == "case", idx] += 2.0
        if gene.startswith("MT-"):
            lam[:, idx] += 0.15
    counts = rng.poisson(lam).astype(np.float32)
    counts[rng.random(counts.shape) < 0.55] = 0
    obs = pd.DataFrame(
        {
            "barcode": [f"CELL{i:03d}-1" for i in range(1, n_cells + 1)],
            "sample_id": sample_id,
            "condition": condition,
        }
    )
    var = pd.DataFrame(
        {
            "gene_id": [f"ENSGDEMO{i:05d}" for i in range(1, n_genes + 1)],
            "gene_symbol": genes,
            "feature_type": "Gene Expression",
            "genome": "GRCh38",
        }
    )
    return sparse.csr_matrix(counts), obs, var


def _write_10x_mtx(counts, obs: pd.DataFrame, var: pd.DataFrame, output_dir: Path) -> None:
    from scipy.io import mmwrite

    output_dir.mkdir(parents=True, exist_ok=True)
    matrix = counts.T.tocoo()
    with gzip.open(output_dir / "matrix.mtx.gz", "wb") as handle:
        mmwrite(handle, matrix)
    with gzip.open(output_dir / "barcodes.tsv.gz", "wt") as handle:
        handle.write("\n".join(obs["barcode"].astype(str)) + "\n")
    with gzip.open(output_dir / "features.tsv.gz", "wt") as handle:
        for row in var.itertuples(index=False):
            handle.write(f"{row.gene_id}\t{row.gene_symbol}\t{row.feature_type}\n")


def _write_10x_h5(counts, obs: pd.DataFrame, var: pd.DataFrame, path: Path) -> None:
    import h5py

    matrix = counts.T.tocsc()
    with h5py.File(path, "w") as h5:
        group = h5.create_group("matrix")
        group.create_dataset("data", data=matrix.data.astype(np.float32), compression="gzip")
        group.create_dataset("indices", data=matrix.indices.astype(np.int64), compression="gzip")
        group.create_dataset("indptr", data=matrix.indptr.astype(np.int64), compression="gzip")
        group.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int64))
        group.create_dataset("barcodes", data=np.asarray(obs["barcode"].astype("S")))
        features = group.create_group("features")
        features.create_dataset("id", data=np.asarray(var["gene_id"].astype("S")))
        features.create_dataset("name", data=np.asarray(var["gene_symbol"].astype("S")))
        features.create_dataset("feature_type", data=np.asarray(var["feature_type"].astype("S")))
        features.create_dataset("genome", data=np.asarray(var["genome"].astype("S")))


INTEGRATION_BACKEND_ID = "scrna.integration.harmony"
INTEGRATION_DEFAULT_BATCH_COLUMN = "batch"


def _resolve_batch_column(config: dict | None) -> str:
    """design.batch_column -> modules.scrna.batch_column -> "batch"."""
    cfg = config if isinstance(config, dict) else {}
    design = cfg.get("design") if isinstance(cfg.get("design"), dict) else {}
    if design.get("batch_column"):
        return str(design["batch_column"])
    modules = cfg.get("modules") if isinstance(cfg.get("modules"), dict) else {}
    scrna_cfg = modules.get("scrna") if isinstance(modules.get("scrna"), dict) else {}
    if scrna_cfg.get("batch_column"):
        return str(scrna_cfg["batch_column"])
    return INTEGRATION_DEFAULT_BATCH_COLUMN


def _resolve_batch_correction_mode(config: dict | None) -> str:
    cfg = config if isinstance(config, dict) else {}
    modules = cfg.get("modules") if isinstance(cfg.get("modules"), dict) else {}
    scrna_cfg = modules.get("scrna") if isinstance(modules.get("scrna"), dict) else {}
    mode = str(scrna_cfg.get("batch_correction") or "none").strip().lower()
    if mode not in {"none", "harmony", "auto"}:
        return "none"
    return mode


def _knn_batch_entropy(embedding: np.ndarray, batches: np.ndarray, k: int = 15) -> float:
    """Mean per-cell batch entropy over its k nearest neighbours (higher = better mixing)."""
    n = embedding.shape[0]
    if n <= 2:
        return 0.0
    k = int(max(1, min(k, n - 1)))
    d2 = ((embedding[:, None, :] - embedding[None, :, :]) ** 2).sum(axis=2)
    np.fill_diagonal(d2, np.inf)
    neighbours = np.argsort(d2, axis=1)[:, :k]
    unique_batches = np.unique(batches)
    entropies = []
    for row in neighbours:
        counts = np.array([(batches[row] == b).sum() for b in unique_batches], dtype=float)
        p = counts[counts > 0] / counts.sum()
        entropies.append(float(-(p * np.log(p)).sum()))
    normaliser = float(np.log(max(2, len(unique_batches))))
    return float(np.mean(entropies) / normaliser) if normaliser > 0 else 0.0


def _knn_ilisi_approx(embedding: np.ndarray, batches: np.ndarray, k: int = 15) -> float:
    """Approximate inverse Simpson index over the neighbourhood batch composition."""
    n = embedding.shape[0]
    if n <= 2:
        return 0.0
    k = int(max(1, min(k, n - 1)))
    d2 = ((embedding[:, None, :] - embedding[None, :, :]) ** 2).sum(axis=2)
    np.fill_diagonal(d2, np.inf)
    neighbours = np.argsort(d2, axis=1)[:, :k]
    unique_batches = np.unique(batches)
    values = []
    for row in neighbours:
        counts = np.array([(batches[row] == b).sum() for b in unique_batches], dtype=float)
        p = counts[counts > 0] / counts.sum()
        values.append(float(1.0 / (p ** 2).sum()))
    return float(np.mean(values))


def _cluster_preservation(before: np.ndarray, after: np.ndarray, seed: int = 7) -> float:
    """Adjusted Rand index between clustering before/after integration."""
    try:
        from sklearn.metrics import adjusted_rand_score
    except ImportError:
        return float("nan")
    return float(adjusted_rand_score(before, after))


def _write_batch_integration_metrics(path: Path, rows: list[dict[str, Any]]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _plot_umap_before_after(before: np.ndarray, after: np.ndarray, batches: np.ndarray, path: Path) -> None:
    import matplotlib.pyplot as plt

    apply_clinical_journal_style()
    unique = list(dict.fromkeys(batches.tolist()))
    cmap = plt.get_cmap("tab10")
    colours = {batch: cmap(idx % 10) for idx, batch in enumerate(unique)}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for ax, coords, subtitle in (
        (axes[0], before, "Before integration"),
        (axes[1], after, "After integration"),
    ):
        for batch in unique:
            mask = batches == batch
            ax.scatter(coords[mask, 0], coords[mask, 1], s=4, color=colours[batch], label=str(batch), linewidth=0)
        ax.set_title(subtitle, fontsize=11)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])
    axes[0].legend(markerscale=3, fontsize=8, title="batch", loc="best")
    fig.tight_layout()
    save_figure(path, close=True)
    plt.close(fig)


def _run_batch_integration(
    *,
    sc: Any,
    adata: Any,
    tables: Path,
    figures: Path,
    config: dict | None,
    random_seed: int,
    analysis_level: str,
) -> tuple[dict[str, Any], dict[str, str]]:
    """Harmony batch integration over X_pca, or an explicit skipped record."""
    backend_id = INTEGRATION_BACKEND_ID
    metrics_path = tables / "batch_integration_metrics.tsv"
    figure_path = figures / "umap_batch_before_after.png"
    status_path = tables / "batch_integration_status.tsv"

    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    batch_column = _resolve_batch_column(config)
    mode = _resolve_batch_correction_mode(config)
    artifacts: dict[str, str] = {
        "batch_integration_metrics": str(metrics_path),
        "batch_integration_status": str(status_path),
        "umap_batch_before_after": str(figure_path),
    }

    if batch_column not in adata.obs:
        return _integration_skip(
            backend_id, analysis_level, f"missing_batch_column:{batch_column}",
            metrics_path, status_path, artifacts,
        )
    batches = adata.obs[batch_column].astype(str).to_numpy()
    n_batches = int(pd.Series(batches).nunique())
    if n_batches < 2:
        return _integration_skip(
            backend_id, analysis_level, f"single_batch_no_integration:{batch_column}",
            metrics_path, status_path, artifacts,
        )
    if mode == "none":
        return _integration_skip(
            backend_id, analysis_level, "batch_correction_disabled:none",
            metrics_path, status_path, artifacts,
        )

    try:
        import harmonypy  # noqa: F401
    except ImportError:
        return _integration_skip(
            backend_id, analysis_level, "dependency_missing:harmonypy",
            metrics_path, status_path, artifacts,
        )

    if "X_pca" not in adata.obsm:
        return _integration_skip(
            backend_id, analysis_level, "missing_pca_embedding",
            metrics_path, status_path, artifacts,
        )

    n_pcs = int(adata.obsm["X_pca"].shape[1])
    # Preserve the pre-integration UMAP/embedding for comparison and reproducibility.
    adata.obsm["X_pca_unintegrated"] = np.asarray(adata.obsm["X_pca"]).copy()
    sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs // 4)), n_pcs=n_pcs)
    sc.tl.umap(adata, random_state=random_seed)
    before_umap = np.asarray(adata.obsm["X_umap"]).copy()
    adata.obsm["X_umap_unintegrated"] = before_umap

    meta = pd.DataFrame({batch_column: batches}, index=adata.obs_names)
    try:
        import harmonypy

        ho = harmonypy.run_harmony(adata.obsm["X_pca"], meta, batch_column, random_state=random_seed, verbose=False)
        corrected = np.asarray(ho.Z_corr)
        if corrected.shape[0] == n_pcs and corrected.shape[1] == adata.n_obs:
            corrected = corrected.T
    except Exception as exc:  # pragma: no cover - depends on optional package internals
        return _integration_skip(
            backend_id, analysis_level, f"harmony_failed:{type(exc).__name__}",
            metrics_path, status_path, artifacts,
        )

    adata.obsm["X_pca_harmony"] = corrected
    adata.uns["ultimate_integration_use_rep"] = "X_pca_harmony"
    cluster_before = _cluster_on_embedding(adata.obsm["X_pca"], random_seed)
    sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs // 4)), use_rep="X_pca_harmony")
    sc.tl.umap(adata, random_state=random_seed)
    after_umap = np.asarray(adata.obsm["X_umap"]).copy()
    cluster_after = _cluster_on_embedding(adata.obsm["X_pca_harmony"], random_seed)

    entropy_before = _knn_batch_entropy(np.asarray(adata.obsm["X_pca_unintegrated"]), batches)
    entropy_after = _knn_batch_entropy(corrected, batches)
    ilisi_before = _knn_ilisi_approx(np.asarray(adata.obsm["X_pca_unintegrated"]), batches)
    ilisi_after = _knn_ilisi_approx(corrected, batches)
    preservation = _cluster_preservation(cluster_before, cluster_after, random_seed)

    _write_batch_integration_metrics(
        metrics_path,
        [
            {
                "backend_id": backend_id,
                "status": "ready",
                "analysis_level": analysis_level,
                "batch_column": batch_column,
                "n_batches": n_batches,
                "n_cells": int(adata.n_obs),
                "n_pcs": n_pcs,
                "batch_entropy_before": round(entropy_before, 6),
                "batch_entropy_after": round(entropy_after, 6),
                "ilisi_approx_before": round(ilisi_before, 6),
                "ilisi_approx_after": round(ilisi_after, 6),
                "cluster_preservation_ari": round(preservation, 6) if preservation == preservation else "",
                "metric_note": "entropy/ilisi are kNN approximations, not scib ilisi",
                "delivery_allowed": False,
            }
        ],
    )
    pd.DataFrame(
        [{"backend_id": backend_id, "status": "ready", "reason": "", "analysis_level": analysis_level,
          "use_rep": "X_pca_harmony", "harmonypy_version": _package_version("harmonypy"), "delivery_allowed": False}]
    ).to_csv(status_path, sep="\t", index=False)
    _plot_umap_before_after(before_umap, after_umap, batches, figure_path)
    adata.obs["batch"] = pd.Categorical(batches)
    return _backend_row(backend_id, "ready", analysis_level, ""), artifacts


def _cluster_on_embedding(embedding: np.ndarray, seed: int) -> np.ndarray:
    try:
        from sklearn.cluster import KMeans

        k = int(min(6, max(2, embedding.shape[0] // 30)))
        return KMeans(n_clusters=k, random_state=seed, n_init=10).fit_predict(np.asarray(embedding))
    except ImportError:  # pragma: no cover
        return np.zeros(embedding.shape[0], dtype=int)


def _package_version(name: str) -> str:
    try:
        module = __import__(name)
    except ImportError:
        return ""
    return str(getattr(module, "__version__", "") or "")


def _integration_skip(
    backend_id: str,
    analysis_level: str,
    reason: str,
    metrics_path: Path,
    status_path: Path,
    artifacts: dict[str, str],
) -> tuple[dict[str, Any], dict[str, str]]:
    """Record an explicit skipped integration, never a silent no-op."""
    _write_batch_integration_metrics(
        metrics_path,
        [
            {
                "backend_id": backend_id,
                "status": "skipped",
                "reason": reason,
                "analysis_level": analysis_level,
                "delivery_allowed": False,
            }
        ],
    )
    pd.DataFrame(
        [{"backend_id": backend_id, "status": "skipped", "reason": reason, "analysis_level": analysis_level, "delivery_allowed": False}]
    ).to_csv(status_path, sep="\t", index=False)
    return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts


def _attach_samplesheet(adata, samplesheet: Path | None) -> None:
    if samplesheet is None or not samplesheet.exists():
        return
    sheet = pd.read_csv(samplesheet, sep=None, engine="python")
    if "barcode" not in sheet.columns:
        return
    sheet = sheet.set_index("barcode")
    for column in sheet.columns:
        mapped = pd.Series(adata.obs_names, index=adata.obs_names).map(sheet[column])
        fallback = adata.obs[column] if column in adata.obs else "unknown"
        adata.obs[column] = mapped.fillna(fallback).to_numpy()


def _qc(sc, adata, figures: Path, tables: Path, analysis_level: str) -> None:
    adata.var["mt"] = pd.Index(adata.var_names.astype(str)).str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    qc = adata.obs[qc_cols].copy()
    qc["analysis_level"] = analysis_level
    qc.to_csv(tables / "qc_metrics.tsv", sep="\t")
    apply_clinical_journal_style()
    import matplotlib.pyplot as plt
    import seaborn as sns

    melted = adata.obs[qc_cols].melt(var_name="metric", value_name="value")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    sns.violinplot(data=melted, x="metric", y="value", cut=0, inner="quartile", ax=ax, color="#7EA6C8")
    ax.set_xlabel("")
    ax.set_ylabel("Value")
    fig.tight_layout()
    save_figure(figures / "qc_violin.png")


def _score_signatures(sc, adata) -> None:
    genes = set(adata.var_names.astype(str))
    rows = []
    for name, gene_set in SIGNATURES.items():
        present = [gene for gene in gene_set if gene in genes]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"{name}_score")
        else:
            adata.obs[f"{name}_score"] = 0.0
        rows.append({"signature": name, "n_input_genes": len(gene_set), "n_present_genes": len(present), "present_genes": ",".join(present)})
    adata.uns["signature_enrichment"] = pd.DataFrame(rows)


def _write_tables(sc, adata, tables: Path, analysis_level: str) -> None:
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    _add_analysis_level(markers, analysis_level).to_csv(tables / "marker_genes.tsv", sep="\t", index=False)
    condition_de_reason = _condition_de_blocker(adata)
    if condition_de_reason:
        _write_condition_de_skip_table(tables / "de_condition.tsv", adata, analysis_level, condition_de_reason)
    else:
        sc.tl.rank_genes_groups(adata, groupby="condition", method="wilcoxon", use_raw=adata.raw is not None)
        condition_de = sc.get.rank_genes_groups_df(adata, group=None)
        _add_analysis_level(condition_de, analysis_level).to_csv(tables / "de_condition.tsv", sep="\t", index=False)
    prop = adata.obs.groupby(["sample_id", "cluster"], observed=False).size().rename("n_cells").reset_index()
    prop["fraction"] = prop["n_cells"] / prop.groupby("sample_id")["n_cells"].transform("sum")
    prop["annotation_status"] = "cluster_only_not_cell_type"
    prop["analysis_level"] = analysis_level
    prop.to_csv(tables / "cell_type_composition.tsv", sep="\t", index=False)
    score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
    if score_cols:
        scores = adata.obs.groupby("cluster", observed=False)[score_cols].mean().reset_index()
        scores["analysis_level"] = analysis_level
        scores.to_csv(tables / "signature_scores_by_cluster.tsv", sep="\t", index=False)
    enrichment = adata.uns["signature_enrichment"].copy()
    enrichment["analysis_level"] = analysis_level
    enrichment.to_csv(tables / "basic_enrichment.tsv", sep="\t", index=False)


def _run_scrublet_backend(
    *,
    counts_adata,
    adata,
    tables: Path,
    figures: Path,
    analysis_level: str,
    random_seed: int,
) -> tuple[dict[str, Any], dict[str, str]]:
    backend_id = "scrna.qc.scrublet"
    score_path = tables / "doublet_scores.tsv"
    summary_path = tables / "doublet_summary.tsv"
    figure_path = figures / "doublet_score_histogram.png"
    n_cells = int(getattr(counts_adata, "n_obs", counts_adata.shape[0]))
    n_features = int(getattr(counts_adata, "n_vars", counts_adata.shape[1]))
    if min(n_cells, n_features) <= 30:
        reason = f"input_too_small_for_scrublet:n_cells={n_cells},n_features={n_features},min_required=31"
        _write_skip_table(score_path, backend_id, analysis_level, reason)
        _write_skip_table(summary_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "skipped", analysis_level, reason), {
            "doublet_scores": str(score_path),
            "doublet_summary": str(summary_path),
        }
    if importlib.util.find_spec("scrublet") is None:
        reason = "dependency_missing:scrublet"
        _write_skip_table(score_path, backend_id, analysis_level, reason)
        _write_skip_table(summary_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "skipped", analysis_level, reason), {
            "doublet_scores": str(score_path),
            "doublet_summary": str(summary_path),
        }
    try:
        import matplotlib.pyplot as plt
        import scrublet as scr

        scores, predicted, n_prin_comps, min_gene_variability_pctl = _scrub_doublets_with_adaptive_pcs(
            scr.Scrublet,
            counts_adata.X,
            random_seed=random_seed,
            n_cells=n_cells,
            n_features=n_features,
        )
        adata.obs["scrublet_score"] = scores
        adata.obs["scrublet_predicted_doublet"] = predicted.astype(bool)
        pd.DataFrame(
            {
                "barcode": adata.obs_names.astype(str),
                "scrublet_score": scores,
                "predicted_doublet": predicted.astype(bool),
                "analysis_level": analysis_level,
            }
        ).to_csv(score_path, sep="\t", index=False)
        pd.DataFrame(
            [
                {
                    "backend_id": backend_id,
                    "status": "ready",
                    "n_cells": int(len(scores)),
                    "predicted_doublets": int(np.asarray(predicted).sum()),
                    "doublet_rate": float(np.asarray(predicted).mean()) if len(scores) else 0.0,
                    "analysis_level": analysis_level,
                    "warning": "doublet calls are model- and threshold-dependent QC flags, not cell type labels",
                    "n_prin_comps": int(n_prin_comps),
                    "min_gene_variability_pctl": float(min_gene_variability_pctl),
                }
            ]
        ).to_csv(summary_path, sep="\t", index=False)
        apply_clinical_journal_style()
        fig, ax = plt.subplots(figsize=(5.8, 4.0))
        ax.hist(scores, bins=30, color="#7EA6C8", edgecolor="white")
        ax.set_xlabel("Scrublet score")
        ax.set_ylabel("Cells")
        fig.tight_layout()
        save_figure(figure_path)
        return _backend_row(backend_id, "ready", analysis_level, ""), {
            "doublet_scores": str(score_path),
            "doublet_summary": str(summary_path),
            "doublet_score_histogram": str(figure_path),
        }
    except Exception as exc:
        if isinstance(exc, ValueError) and "n_components" in str(exc):
            reason = f"input_too_small_for_scrublet_after_filtering:{exc}"
            _write_skip_table(score_path, backend_id, analysis_level, reason)
            _write_skip_table(summary_path, backend_id, analysis_level, reason)
            return _backend_row(backend_id, "skipped", analysis_level, reason), {
                "doublet_scores": str(score_path),
                "doublet_summary": str(summary_path),
            }
        reason = f"backend_failed:{type(exc).__name__}:{exc}"
        _write_skip_table(score_path, backend_id, analysis_level, reason)
        _write_skip_table(summary_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "failed", analysis_level, reason), {
            "doublet_scores": str(score_path),
            "doublet_summary": str(summary_path),
        }


def _scrub_doublets_with_adaptive_pcs(scrublet_cls, counts_matrix, *, random_seed: int, n_cells: int, n_features: int):
    max_comps = max(2, min(30, n_cells - 1, n_features - 1))
    pc_candidates = []
    for value in (max_comps, 20, 15, 10, 5, 2):
        value = max(2, min(value, max_comps))
        if value not in pc_candidates:
            pc_candidates.append(value)
    variability_candidates = (85, 75, 50, 25, 10, 5, 0)
    last_error: Exception | None = None
    for min_gene_variability_pctl in variability_candidates:
        for n_prin_comps in pc_candidates:
            try:
                scrub = scrublet_cls(counts_matrix, random_state=random_seed)
                scores, predicted = scrub.scrub_doublets(
                    n_prin_comps=n_prin_comps,
                    min_gene_variability_pctl=min_gene_variability_pctl,
                    verbose=False,
                )
                return scores, predicted, n_prin_comps, min_gene_variability_pctl
            except (IndexError, ValueError) as exc:
                if not _is_retryable_scrublet_parameter_error(exc):
                    raise
                last_error = exc
    raise last_error or ValueError("scrublet failed to select a valid n_prin_comps")


def _is_retryable_scrublet_parameter_error(exc: Exception) -> bool:
    message = str(exc)
    return (
        "n_components" in message
        or "0 feature(s)" in message
        or ("index -1 is out of bounds" in message and "size 0" in message)
    )


def _run_celltypist_backend(
    *,
    adata,
    tables: Path,
    analysis_level: str,
    model_path: Path | None,
    figures: Path | None = None,
    objects: Path | None = None,
    manifests: Path | None = None,
    reports: Path | None = None,
    majority_voting: bool = True,
    cluster_key: str = "cluster",
    uncertainty_threshold: float = 0.5,
) -> tuple[dict[str, Any], dict[str, str]]:
    backend_id = "scrna.annotation.celltypist"
    annotation_path = tables / "celltypist_annotation.tsv"
    confidence_path = tables / "annotation_confidence.tsv"
    decision_path = tables / "celltypist_decision_scores.tsv"
    feature_summary_path = tables / "celltypist_feature_match_summary.tsv"
    missing_model_path = tables / "celltypist_model_features_missing_from_input.tsv"
    unmatched_input_path = tables / "celltypist_input_genes_not_in_model.tsv"
    duplicate_gene_path = tables / "celltypist_duplicate_genes.tsv"
    label_distribution_path = tables / "celltypist_label_distribution.tsv"
    cluster_crosstab_path = tables / "celltypist_cluster_crosstab.tsv"
    marker_review_path = tables / "celltypist_marker_review.tsv"
    uncertainty_path = tables / "celltypist_uncertain_cells.tsv"
    warning_path = tables / "annotation_warning.tsv"
    manifest_path = (manifests or tables) / "celltypist_annotation_manifest.json"
    annotated_h5ad_path = (objects / "scrna_celltypist_annotated.h5ad") if objects is not None else None
    umap_label_path = (figures / "celltypist_umap_label.png") if figures is not None else None
    umap_score_path = (figures / "celltypist_umap_max_probability.png") if figures is not None else None
    marker_dotplot_path = (figures / "celltypist_marker_dotplot.png") if figures is not None else None
    report_path = (reports / "celltypist_annotation_report.md") if reports is not None else None
    if importlib.util.find_spec("celltypist") is None:
        reason = "dependency_missing:celltypist"
        _write_skip_table(annotation_path, backend_id, analysis_level, reason)
        _write_skip_table(confidence_path, backend_id, analysis_level, reason)
        _write_annotation_warning_table(warning_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "skipped", analysis_level, reason), _celltypist_artifacts(annotation_path, confidence_path, warning_path)
    if model_path is None or not model_path.exists():
        reason = "model_not_configured:celltypist_model_required_to_avoid_implicit_download"
        _write_skip_table(annotation_path, backend_id, analysis_level, reason)
        _write_skip_table(confidence_path, backend_id, analysis_level, reason)
        _write_annotation_warning_table(warning_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "skipped", analysis_level, reason), _celltypist_artifacts(annotation_path, confidence_path, warning_path)
    try:
        import celltypist
        from celltypist import models as celltypist_models

        annotation_adata = _full_gene_expression_adata(adata).copy()
        _validate_celltypist_expression_input(annotation_adata)
        model = celltypist_models.Model.load(str(model_path))
        model_features = _celltypist_model_features(model)
        feature_report = _write_celltypist_feature_match_tables(
            annotation_adata=annotation_adata,
            model_features=model_features,
            model_path=model_path,
            feature_summary_path=feature_summary_path,
            missing_model_path=missing_model_path,
            unmatched_input_path=unmatched_input_path,
            duplicate_gene_path=duplicate_gene_path,
            analysis_level=analysis_level,
        )
        if feature_report["n_matched_model_features"] == 0:
            raise ValueError("no_matching_model_genes_between_input_and_celltypist_model")
        voting_key = cluster_key if majority_voting and cluster_key in annotation_adata.obs.columns else None
        result = celltypist.annotate(
            annotation_adata,
            model=model,
            majority_voting=bool(majority_voting),
            over_clustering=voting_key,
        )
        predictions_indexed = result.predicted_labels.copy()
        probability_matrix = getattr(result, "probability_matrix", None)
        decision_matrix = getattr(result, "decision_matrix", None)
        annotation_table = _celltypist_annotation_table(
            predictions_indexed=predictions_indexed,
            probability_matrix=probability_matrix,
            decision_matrix=decision_matrix,
            cell_index=annotation_adata.obs_names,
            obs=annotation_adata.obs,
            cluster_key=cluster_key,
            uncertainty_threshold=uncertainty_threshold,
            analysis_level=analysis_level,
        )
        annotation_table.to_csv(annotation_path, sep="\t", index=False)
        _write_celltypist_matrix(probability_matrix, confidence_path, backend_id, analysis_level, "probability_matrix_not_available")
        _write_celltypist_matrix(decision_matrix, decision_path, backend_id, analysis_level, "decision_matrix_not_available")
        _write_celltypist_label_distribution(annotation_table, label_distribution_path, analysis_level)
        _write_celltypist_cluster_crosstab(annotation_table, cluster_crosstab_path, analysis_level)
        annotation_table.loc[annotation_table["uncertain"].astype(bool)].to_csv(uncertainty_path, sep="\t", index=False)
        marker_review = _write_celltypist_marker_review(
            adata=annotation_adata,
            annotation_table=annotation_table,
            marker_review_path=marker_review_path,
            analysis_level=analysis_level,
        )
        _write_celltypist_figures(
            adata=annotation_adata,
            annotation_table=annotation_table,
            umap_label_path=umap_label_path,
            umap_score_path=umap_score_path,
            marker_dotplot_path=marker_dotplot_path,
            marker_review=marker_review,
        )
        has_majority_labels = bool(annotation_table["celltypist_majority_label"].astype(str).ne("").any())
        final_label_column = "celltypist_majority_label" if has_majority_labels else "celltypist_raw_label"
        for column in [
            "celltypist_raw_label",
            "celltypist_majority_label",
            "celltypist_label",
            "celltypist_max_probability",
            "celltypist_raw_label_probability",
            "celltypist_uncertain",
        ]:
            if column in annotation_table.columns:
                values = annotation_table.set_index("barcode")[column].reindex(annotation_adata.obs_names)
                if column.endswith("_label") or column == "celltypist_label":
                    annotation_adata.obs[column] = pd.Categorical(values.fillna("unassigned").astype(str).to_numpy())
                else:
                    annotation_adata.obs[column] = values.to_numpy()
        annotation_adata.obs["celltypist_label_source"] = f"celltypist:{model_path.name}"
        annotation_adata.uns["ultimate_celltypist"] = {
            **feature_report,
            "model_path": str(model_path),
            "model_sha256": _sha256_file(model_path),
            "celltypist_version": getattr(celltypist, "__version__", "unknown"),
            "majority_voting_requested": bool(majority_voting),
            "majority_voting_used": bool(voting_key),
            "over_clustering": voting_key or "",
            "uncertainty_threshold": float(uncertainty_threshold),
            "uncertainty_rule": "exploratory:max_probability_below_threshold; not a calibrated accuracy estimate",
            "final_label_column": final_label_column,
        }
        if annotated_h5ad_path is not None:
            annotation_adata.write_h5ad(annotated_h5ad_path)
        manifest = _write_celltypist_annotation_manifest(
            manifest_path=manifest_path,
            annotation_adata=annotation_adata,
            annotation_table=annotation_table,
            feature_report=feature_report,
            model_path=model_path,
            analysis_level=analysis_level,
            majority_voting=majority_voting,
            voting_key=voting_key,
            uncertainty_threshold=uncertainty_threshold,
            artifacts=_celltypist_artifacts(
                annotation_path,
                confidence_path,
                warning_path,
                decision_path=decision_path,
                feature_summary_path=feature_summary_path,
                missing_model_path=missing_model_path,
                unmatched_input_path=unmatched_input_path,
                duplicate_gene_path=duplicate_gene_path,
                label_distribution_path=label_distribution_path,
                cluster_crosstab_path=cluster_crosstab_path,
                marker_review_path=marker_review_path,
                uncertainty_path=uncertainty_path,
                annotated_h5ad_path=annotated_h5ad_path,
                manifest_path=manifest_path,
                report_path=report_path,
                umap_label_path=umap_label_path,
                umap_score_path=umap_score_path,
                marker_dotplot_path=marker_dotplot_path,
            ),
        )
        if report_path is not None:
            _write_celltypist_report(report_path, manifest)
        _write_annotation_warning_table(
            warning_path,
            backend_id,
            analysis_level,
            "CellTypist labels depend on the supplied model/reference; max probability is a model score, not a calibrated probability that the label is correct.",
        )
        return _backend_row(backend_id, "ready", analysis_level, ""), _celltypist_artifacts(
            annotation_path,
            confidence_path,
            warning_path,
            decision_path=decision_path,
            feature_summary_path=feature_summary_path,
            missing_model_path=missing_model_path,
            unmatched_input_path=unmatched_input_path,
            duplicate_gene_path=duplicate_gene_path,
            label_distribution_path=label_distribution_path,
            cluster_crosstab_path=cluster_crosstab_path,
            marker_review_path=marker_review_path,
            uncertainty_path=uncertainty_path,
            annotated_h5ad_path=annotated_h5ad_path,
            manifest_path=manifest_path,
            report_path=report_path,
            umap_label_path=umap_label_path,
            umap_score_path=umap_score_path,
            marker_dotplot_path=marker_dotplot_path,
        )
    except Exception as exc:
        reason = f"backend_failed:{type(exc).__name__}:{exc}"
        _write_skip_table(annotation_path, backend_id, analysis_level, reason)
        _write_skip_table(confidence_path, backend_id, analysis_level, reason)
        _write_annotation_warning_table(warning_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "failed", analysis_level, reason), _celltypist_artifacts(annotation_path, confidence_path, warning_path)


def _validate_celltypist_expression_input(adata) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("empty_celltypist_input")
    if not adata.obs_names.is_unique:
        raise ValueError("celltypist_input_obs_names_not_unique")
    if not adata.var_names.is_unique:
        raise ValueError("celltypist_input_var_names_not_unique")
    if _has_negative_values(adata.X):
        raise ValueError("celltypist_input_contains_negative_values_expected_log_normalized_expression")
    from scipy import sparse

    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if values.size and not np.isfinite(values).all():
        raise ValueError("celltypist_input_contains_non_finite_values")


def _celltypist_model_features(model) -> list[str]:
    features = getattr(model, "features", None)
    if features is None:
        classifier = getattr(model, "classifier", None)
        features = getattr(classifier, "feature_names_in_", None)
    if features is None:
        raise ValueError("celltypist_model_features_not_available")
    return [str(gene) for gene in list(features)]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_celltypist_feature_match_tables(
    *,
    annotation_adata,
    model_features: list[str],
    model_path: Path,
    feature_summary_path: Path,
    missing_model_path: Path,
    unmatched_input_path: Path,
    duplicate_gene_path: Path,
    analysis_level: str,
) -> dict[str, Any]:
    input_genes = pd.Index(annotation_adata.var_names.astype(str))
    input_gene_set = set(input_genes)
    model_gene_set = set(model_features)
    matched_model_features = [gene for gene in model_features if gene in input_gene_set]
    missing_model_features = [gene for gene in model_features if gene not in input_gene_set]
    unmatched_input_genes = [gene for gene in input_genes if gene not in model_gene_set]
    duplicate_input_genes = input_genes[input_genes.duplicated()].tolist()
    duplicate_model_features = pd.Index(model_features)[pd.Index(model_features).duplicated()].tolist()
    summary = {
        "backend_id": "scrna.annotation.celltypist",
        "model_path": str(model_path),
        "model_file": model_path.name,
        "model_sha256": _sha256_file(model_path),
        "n_input_cells": int(annotation_adata.n_obs),
        "n_input_genes": int(annotation_adata.n_vars),
        "n_model_features": int(len(model_features)),
        "n_matched_model_features": int(len(matched_model_features)),
        "model_feature_match_fraction": float(len(matched_model_features) / len(model_features)) if model_features else 0.0,
        "n_missing_model_features": int(len(missing_model_features)),
        "n_input_genes_not_in_model": int(len(unmatched_input_genes)),
        "n_duplicate_input_gene_names": int(len(duplicate_input_genes)),
        "n_duplicate_model_features": int(len(duplicate_model_features)),
        "gene_identifier_basis": "AnnData.var_names matched directly to CellTypist model features; no case conversion or ID mapping was applied.",
        "input_matrix_semantics": _celltypist_input_semantics(annotation_adata),
        "analysis_level": analysis_level,
    }
    pd.DataFrame([summary]).to_csv(feature_summary_path, sep="\t", index=False)
    pd.DataFrame({"gene": missing_model_features, "analysis_level": analysis_level}).to_csv(missing_model_path, sep="\t", index=False)
    pd.DataFrame({"gene": unmatched_input_genes, "analysis_level": analysis_level}).to_csv(unmatched_input_path, sep="\t", index=False)
    pd.DataFrame(
        [
            {"source": "input", "gene": gene, "analysis_level": analysis_level}
            for gene in duplicate_input_genes
        ]
        + [
            {"source": "model", "gene": gene, "analysis_level": analysis_level}
            for gene in duplicate_model_features
        ]
    ).to_csv(duplicate_gene_path, sep="\t", index=False)
    return summary


def _celltypist_input_semantics(adata) -> str:
    semantics = getattr(adata, "uns", {}).get("ultimate_matrix_semantics")
    if isinstance(semantics, dict) and semantics.get("X"):
        return str(semantics["X"])
    if _has_log1p_metadata(adata):
        return "AnnData.X has log1p metadata; expected normalized/log expression"
    return "AnnData.X passed through as provided; no normalization or log transform was applied by this annotation entrypoint"


def _celltypist_annotation_table(
    *,
    predictions_indexed: pd.DataFrame,
    probability_matrix: pd.DataFrame | None,
    decision_matrix: pd.DataFrame | None,
    cell_index: pd.Index,
    obs: pd.DataFrame,
    cluster_key: str,
    uncertainty_threshold: float,
    analysis_level: str,
) -> pd.DataFrame:
    cell_ids = pd.Index(cell_index.astype(str), name="barcode")
    predictions = predictions_indexed.copy()
    predictions.index = predictions.index.astype(str)
    if not predictions.index.is_unique:
        raise ValueError("celltypist_prediction_index_not_unique")
    missing = cell_ids.difference(predictions.index)
    extra = predictions.index.difference(cell_ids)
    if len(missing) or len(extra):
        raise ValueError(f"celltypist_prediction_cell_ids_not_aligned:missing={len(missing)},extra={len(extra)}")
    predictions = predictions.loc[cell_ids]
    raw_labels = predictions["predicted_labels"].astype(str) if "predicted_labels" in predictions.columns else predictions.iloc[:, 0].astype(str)
    majority_labels = (
        predictions["majority_voting"].astype(str)
        if "majority_voting" in predictions.columns
        else pd.Series("", index=cell_ids, dtype="object")
    )
    final_labels = majority_labels.where(majority_labels != "", raw_labels)
    probabilities = _celltypist_aligned_matrix(probability_matrix, cell_ids)
    decisions = _celltypist_aligned_matrix(decision_matrix, cell_ids)
    if probabilities is not None and not probabilities.empty:
        max_probability = probabilities.max(axis=1)
        max_probability_label = probabilities.idxmax(axis=1).astype(str)
        raw_label_probability = _lookup_celltypist_scores(probabilities, raw_labels)
        majority_label_probability = _lookup_celltypist_scores(probabilities, majority_labels.replace("", np.nan))
    else:
        max_probability = pd.Series(np.nan, index=cell_ids)
        max_probability_label = pd.Series("", index=cell_ids, dtype="object")
        raw_label_probability = pd.Series(np.nan, index=cell_ids)
        majority_label_probability = pd.Series(np.nan, index=cell_ids)
    if decisions is not None and not decisions.empty:
        raw_label_decision_score = _lookup_celltypist_scores(decisions, raw_labels)
        majority_label_decision_score = _lookup_celltypist_scores(decisions, majority_labels.replace("", np.nan))
    else:
        raw_label_decision_score = pd.Series(np.nan, index=cell_ids)
        majority_label_decision_score = pd.Series(np.nan, index=cell_ids)
    table = pd.DataFrame(
        {
            "barcode": cell_ids,
            "celltypist_raw_label": raw_labels.to_numpy(),
            "celltypist_majority_label": majority_labels.to_numpy(),
            "celltypist_label": final_labels.to_numpy(),
            "celltypist_max_probability": max_probability.to_numpy(),
            "celltypist_max_probability_label": max_probability_label.to_numpy(),
            "celltypist_raw_label_probability": raw_label_probability.to_numpy(),
            "celltypist_majority_label_probability": majority_label_probability.to_numpy(),
            "celltypist_raw_label_decision_score": raw_label_decision_score.to_numpy(),
            "celltypist_majority_label_decision_score": majority_label_decision_score.to_numpy(),
        }
    )
    if cluster_key in obs.columns:
        table[cluster_key] = obs.loc[cell_ids, cluster_key].astype(str).to_numpy()
    table["uncertainty_threshold"] = float(uncertainty_threshold)
    table["uncertain"] = table["celltypist_max_probability"].isna() | (table["celltypist_max_probability"] < float(uncertainty_threshold))
    table["uncertainty_rule"] = "exploratory:max_probability_below_threshold"
    table["annotation_status"] = np.where(table["uncertain"], "celltypist_prediction_uncertain", "celltypist_model_prediction")
    table["analysis_level"] = analysis_level
    return table


def _celltypist_aligned_matrix(matrix: pd.DataFrame | None, cell_ids: pd.Index) -> pd.DataFrame | None:
    if matrix is None:
        return None
    frame = matrix.copy()
    frame.index = frame.index.astype(str)
    missing = cell_ids.difference(frame.index)
    extra = frame.index.difference(cell_ids)
    if len(missing) or len(extra):
        raise ValueError(f"celltypist_matrix_cell_ids_not_aligned:missing={len(missing)},extra={len(extra)}")
    return frame.loc[cell_ids]


def _lookup_celltypist_scores(matrix: pd.DataFrame, labels: pd.Series) -> pd.Series:
    values: list[float] = []
    for barcode, label in labels.items():
        if pd.isna(label) or str(label) not in matrix.columns:
            values.append(np.nan)
        else:
            values.append(float(matrix.at[barcode, str(label)]))
    return pd.Series(values, index=labels.index)


def _write_celltypist_matrix(
    matrix: pd.DataFrame | None,
    path: Path,
    backend_id: str,
    analysis_level: str,
    missing_reason: str,
) -> None:
    if matrix is None:
        _write_skip_table(path, backend_id, analysis_level, missing_reason)
        return
    frame = matrix.copy()
    frame.index.name = "barcode"
    frame.reset_index().to_csv(path, sep="\t", index=False)


def _write_celltypist_label_distribution(annotation_table: pd.DataFrame, path: Path, analysis_level: str) -> None:
    rows = []
    for column in ("celltypist_raw_label", "celltypist_majority_label", "celltypist_label"):
        if column not in annotation_table.columns:
            continue
        counts = annotation_table[column].replace("", "not_available").value_counts(dropna=False)
        for label, count in counts.items():
            rows.append({"label_source": column, "label": str(label), "n_cells": int(count), "analysis_level": analysis_level})
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _write_celltypist_cluster_crosstab(annotation_table: pd.DataFrame, path: Path, analysis_level: str) -> None:
    if "cluster" not in annotation_table.columns:
        pd.DataFrame(
            [{"status": "skipped", "reason": "cluster_column_not_available", "analysis_level": analysis_level}]
        ).to_csv(path, sep="\t", index=False)
        return
    crosstab = pd.crosstab(annotation_table["cluster"].astype(str), annotation_table["celltypist_label"].astype(str))
    crosstab.index.name = "cluster"
    crosstab.reset_index().assign(analysis_level=analysis_level).to_csv(path, sep="\t", index=False)


def _write_celltypist_marker_review(
    *,
    adata,
    annotation_table: pd.DataFrame,
    marker_review_path: Path,
    analysis_level: str,
) -> pd.DataFrame:
    from scipy import sparse

    label_series = annotation_table.set_index("barcode")["celltypist_label"].astype(str).reindex(adata.obs_names.astype(str))
    genes = pd.Index(adata.var_names.astype(str))
    rows = []
    for marker_panel, marker_genes in CELLTYPIST_MARKER_PANELS.items():
        present = [gene for gene in marker_genes if gene in genes]
        for label in sorted(label_series.dropna().unique()):
            positions = np.where(label_series.to_numpy() == label)[0]
            for gene in present:
                column = genes.get_loc(gene)
                vector = adata.X[positions, column]
                values = vector.toarray().ravel() if sparse.issparse(vector) else np.asarray(vector).ravel()
                rows.append(
                    {
                        "celltypist_label": label,
                        "marker_panel": marker_panel,
                        "marker_gene": gene,
                        "n_cells": int(len(positions)),
                        "mean_log_expression": float(np.nanmean(values)) if values.size else 0.0,
                        "fraction_detected": float(np.mean(values > 0)) if values.size else 0.0,
                        "analysis_level": analysis_level,
                    }
                )
        if not present:
            rows.append(
                {
                    "celltypist_label": "all",
                    "marker_panel": marker_panel,
                    "marker_gene": "",
                    "n_cells": int(adata.n_obs),
                    "mean_log_expression": np.nan,
                    "fraction_detected": np.nan,
                    "analysis_level": analysis_level,
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(marker_review_path, sep="\t", index=False)
    return frame


def _write_celltypist_figures(
    *,
    adata,
    annotation_table: pd.DataFrame,
    umap_label_path: Path | None,
    umap_score_path: Path | None,
    marker_dotplot_path: Path | None,
    marker_review: pd.DataFrame,
) -> None:
    if umap_label_path is None and umap_score_path is None and marker_dotplot_path is None:
        return
    import matplotlib.pyplot as plt
    import seaborn as sns

    apply_clinical_journal_style()
    if "X_umap" in adata.obsm and umap_label_path is not None:
        coords = np.asarray(adata.obsm["X_umap"])
        labels = annotation_table.set_index("barcode")["celltypist_label"].astype(str).reindex(adata.obs_names.astype(str)).fillna("unassigned")
        codes, uniques = pd.factorize(labels)
        fig, ax = plt.subplots(figsize=(6.2, 5.2))
        ax.scatter(coords[:, 0], coords[:, 1], c=codes, s=6, cmap="tab20", linewidths=0)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title(f"CellTypist labels ({len(uniques)} classes)")
        fig.tight_layout()
        save_figure(umap_label_path)
    if "X_umap" in adata.obsm and umap_score_path is not None:
        coords = np.asarray(adata.obsm["X_umap"])
        scores = annotation_table.set_index("barcode")["celltypist_max_probability"].reindex(adata.obs_names.astype(str)).astype(float)
        fig, ax = plt.subplots(figsize=(6.2, 5.2))
        points = ax.scatter(coords[:, 0], coords[:, 1], c=scores, s=6, cmap="viridis", linewidths=0)
        fig.colorbar(points, ax=ax, label="Max CellTypist probability")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title("CellTypist score distribution")
        fig.tight_layout()
        save_figure(umap_score_path)
    if marker_dotplot_path is not None and not marker_review.empty:
        label_counts = annotation_table["celltypist_label"].astype(str).value_counts().head(20).index
        summary = (
            marker_review[marker_review["celltypist_label"].isin(label_counts)]
            .groupby(["celltypist_label", "marker_panel"], observed=False)["mean_log_expression"]
            .mean()
            .reset_index()
        )
        if not summary.empty:
            matrix = summary.pivot(index="celltypist_label", columns="marker_panel", values="mean_log_expression").fillna(0.0)
            fig_height = max(3.0, 0.28 * len(matrix) + 1.4)
            fig, ax = plt.subplots(figsize=(7.2, fig_height))
            sns.heatmap(matrix, cmap="mako", ax=ax, cbar_kws={"label": "Mean log expression"})
            ax.set_xlabel("Marker panel")
            ax.set_ylabel("CellTypist label")
            fig.tight_layout()
            save_figure(marker_dotplot_path)


def _write_celltypist_annotation_manifest(
    *,
    manifest_path: Path,
    annotation_adata,
    annotation_table: pd.DataFrame,
    feature_report: dict[str, Any],
    model_path: Path,
    analysis_level: str,
    majority_voting: bool,
    voting_key: str | None,
    uncertainty_threshold: float,
    artifacts: dict[str, str],
) -> dict[str, Any]:
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "backend_id": "scrna.annotation.celltypist",
        "status": "ready",
        "analysis_level": analysis_level,
        "n_cells": int(annotation_adata.n_obs),
        "n_genes": int(annotation_adata.n_vars),
        "model_path": str(model_path),
        "model_sha256": _sha256_file(model_path),
        "feature_match": feature_report,
        "raw_label_count": int(annotation_table["celltypist_raw_label"].nunique()),
        "final_label_count": int(annotation_table["celltypist_label"].nunique()),
        "uncertain_cells": int(annotation_table["uncertain"].sum()),
        "uncertainty_threshold": float(uncertainty_threshold),
        "uncertainty_rule": "exploratory:max_probability_below_threshold; max probability is not calibrated label correctness.",
        "majority_voting_requested": bool(majority_voting),
        "majority_voting_used": bool(voting_key),
        "over_clustering": voting_key or "",
        "input_contract": "AnnData.X must be normalized/log expression with gene symbols in var_names; layers['counts'] is preserved but not used for CellTypist prediction.",
        "artifacts": artifacts,
        "slurm": {
            "job_id": os.environ.get("SLURM_JOB_ID", ""),
            "job_name": os.environ.get("SLURM_JOB_NAME", ""),
            "submit_dir": os.environ.get("SLURM_SUBMIT_DIR", ""),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def _write_celltypist_report(path: Path, manifest: dict[str, Any]) -> None:
    lines = [
        "# CellTypist annotation validation",
        "",
        f"- status: `{manifest['status']}`",
        f"- analysis_level: `{manifest['analysis_level']}`",
        f"- model: `{manifest['model_path']}`",
        f"- model_sha256: `{manifest['model_sha256']}`",
        f"- input cells × genes: `{manifest['n_cells']} × {manifest['n_genes']}`",
        f"- model features matched: `{manifest['feature_match']['n_matched_model_features']} / {manifest['feature_match']['n_model_features']}`",
        f"- raw labels: `{manifest['raw_label_count']}`",
        f"- final labels: `{manifest['final_label_count']}`",
        f"- uncertain cells: `{manifest['uncertain_cells']}` at exploratory threshold `{manifest['uncertainty_threshold']}`",
        f"- majority voting used: `{manifest['majority_voting_used']}` via `{manifest['over_clustering'] or 'not_used'}`",
        "",
        "CellTypist score values are model outputs used for review and filtering; this report does not treat them as calibrated probabilities that labels are correct.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _attach_celltypist_embedding(adata, embedding):
    adata = adata.copy()
    missing = pd.Index(adata.obs_names.astype(str)).difference(pd.Index(embedding.obs_names.astype(str)))
    if len(missing):
        raise ValueError(f"embedding_h5ad_missing_input_cells:{len(missing)}")
    embedding = embedding[adata.obs_names].copy()
    for key in ("X_umap", "X_pca"):
        if key in embedding.obsm:
            adata.obsm[key] = embedding.obsm[key].copy()
    for column in ("cluster", "leiden", "sample_id", "condition"):
        if column not in adata.obs.columns and column in embedding.obs.columns:
            adata.obs[column] = embedding.obs[column].reindex(adata.obs_names).to_numpy()
    return adata


def _run_functional_backend(*, adata, tables: Path, figures: Path, analysis_level: str) -> tuple[dict[str, Any], dict[str, str]]:
    backend_id = "scrna.functional.decoupler_gseapy"
    signature_path = tables / "signature_scores.tsv"
    pathway_path = tables / "pathway_activity.tsv"
    tf_path = tables / "tf_activity.tsv"
    status_path = tables / "functional_backend_status.tsv"
    figure_path = figures / "signature_heatmap.png"
    artifacts = {
        "signature_scores": str(signature_path),
        "pathway_activity": str(pathway_path),
        "tf_activity": str(tf_path),
        "functional_backend_status": str(status_path),
        "signature_heatmap": str(figure_path),
    }
    missing = [package for package in ("decoupler", "gseapy") if importlib.util.find_spec(package) is None]
    if missing:
        reason = "dependency_missing:" + ",".join(missing)
        for path in (signature_path, pathway_path, tf_path):
            _write_skip_table(path, backend_id, analysis_level, reason)
        _write_functional_status(status_path, backend_id, analysis_level, "skipped", reason)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    try:
        import decoupler as dc
        import gseapy as gp
        import importlib.metadata as md

        cluster_expr = _cluster_expression_matrix(adata)
        signature_sets = _present_gene_sets(cluster_expr.columns, SIGNATURES, min_genes=2)
        tf_sets = _present_gene_sets(cluster_expr.columns, TF_SIGNATURES, min_genes=2)
        if not signature_sets:
            reason = "insufficient_gene_overlap:signature_sets"
            for path in (signature_path, pathway_path, tf_path):
                _write_skip_table(path, backend_id, analysis_level, reason)
            _write_functional_status(status_path, backend_id, analysis_level, "skipped", reason)
            return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts

        _write_signature_scores_table(tables, signature_path, analysis_level, backend_id)
        pathway_rows = []
        ssgsea = gp.ssgsea(
            data=cluster_expr.T,
            gene_sets=signature_sets,
            min_size=2,
            max_size=max(50, int(cluster_expr.shape[1])),
            outdir=None,
            no_plot=True,
            threads=1,
            seed=123,
            verbose=False,
        )
        gseapy_scores = ssgsea.res2d.copy()
        gseapy_scores["backend_id"] = backend_id
        gseapy_scores["method"] = "gseapy_ssgsea"
        gseapy_scores["analysis_level"] = analysis_level
        pathway_rows.append(gseapy_scores)
        decoupler_net = _gene_set_network(signature_sets)
        decoupler_scores = dc.mt.aucell(cluster_expr, decoupler_net, tmin=2, raw=False, verbose=False)
        pathway_rows.append(
            _activity_result_to_long(
                decoupler_scores,
                backend_id=backend_id,
                method="decoupler_aucell",
                value_name="activity_score",
                analysis_level=analysis_level,
            )
        )
        pd.concat(pathway_rows, ignore_index=True, sort=False).to_csv(pathway_path, sep="\t", index=False)

        if tf_sets:
            tf_scores = dc.mt.ulm(cluster_expr, _gene_set_network(tf_sets), tmin=2, raw=False, verbose=False)
            tf_activity = _activity_result_to_long(
                tf_scores,
                backend_id=backend_id,
                method="decoupler_ulm_curated_tf_like_signatures",
                value_name="tf_activity_score",
                analysis_level=analysis_level,
            )
            tf_activity["warning"] = "TF activity uses curated expression-derived target signatures; it is not experimental TF binding evidence."
            tf_activity.to_csv(tf_path, sep="\t", index=False)
        else:
            _write_skip_table(tf_path, backend_id, analysis_level, "insufficient_gene_overlap:tf_signature_sets")

        _write_functional_status(
            status_path,
            backend_id,
            analysis_level,
            "ready",
            "",
            decoupler_version=md.version("decoupler"),
            gseapy_version=md.version("gseapy"),
        )
        return _backend_row(backend_id, "ready", analysis_level, ""), artifacts
    except Exception as exc:
        reason = f"backend_failed:{type(exc).__name__}:{exc}"
        for path in (signature_path, pathway_path, tf_path):
            _write_skip_table(path, backend_id, analysis_level, reason)
        _write_functional_status(status_path, backend_id, analysis_level, "failed", reason)
        return _backend_row(backend_id, "failed", analysis_level, reason), artifacts


def _run_liana_backend(
    *,
    adata,
    tables: Path,
    figures: Path,
    analysis_level: str,
    random_seed: int,
) -> tuple[dict[str, Any], dict[str, str]]:
    backend_id = "scrna.communication.liana"
    interactions_path = tables / "liana_interactions.tsv"
    network_path = tables / "communication_network.tsv"
    status_path = tables / "communication_backend_status.tsv"
    figure_path = figures / "communication_dotplot.png"
    artifacts = {
        "liana_interactions": str(interactions_path),
        "communication_network": str(network_path),
        "communication_backend_status": str(status_path),
        "communication_dotplot": str(figure_path),
    }
    if importlib.util.find_spec("liana") is None:
        reason = "dependency_missing:liana"
        _write_liana_skip_outputs(interactions_path, network_path, status_path, figure_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    groupby, reason = _communication_groupby(adata)
    if reason:
        _write_liana_skip_outputs(interactions_path, network_path, status_path, figure_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    try:
        import importlib.metadata as md
        import liana as li

        work = adata.copy()
        li.mt.rank_aggregate(
            work,
            groupby=groupby,
            resource_name="consensus",
            expr_prop=0.2,
            min_cells=10,
            n_perms=20,
            n_jobs=1,
            seed=random_seed,
            use_raw=True,
            key_added="liana_res",
            inplace=True,
            verbose=False,
        )
        interactions = pd.DataFrame(work.uns.get("liana_res", pd.DataFrame())).copy()
        if interactions.empty:
            reason = "empty_result:liana_res"
            _write_liana_skip_outputs(interactions_path, network_path, status_path, figure_path, backend_id, analysis_level, reason)
            return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
        interactions["backend_id"] = backend_id
        interactions["method"] = "liana_rank_aggregate_consensus"
        interactions["groupby"] = groupby
        interactions["analysis_level"] = analysis_level
        interactions["warning"] = "Ligand-receptor communication is an expression-derived statistical inference, not direct mechanism proof."
        interactions.to_csv(interactions_path, sep="\t", index=False)
        network = _liana_network_summary(interactions, backend_id, groupby, analysis_level)
        network.to_csv(network_path, sep="\t", index=False)
        _write_liana_dotplot(interactions, figure_path)
        _write_liana_status(
            status_path,
            backend_id,
            analysis_level,
            "ready",
            "",
            groupby=groupby,
            liana_version=md.version("liana"),
            n_interactions=int(interactions.shape[0]),
            n_cell_groups=int(work.obs[groupby].astype(str).nunique()),
        )
        return _backend_row(backend_id, "ready", analysis_level, ""), artifacts
    except Exception as exc:
        reason = f"backend_failed:{type(exc).__name__}:{exc}"
        _write_liana_skip_outputs(interactions_path, network_path, status_path, figure_path, backend_id, analysis_level, reason)
        return _backend_row(backend_id, "failed", analysis_level, reason), artifacts


def _run_cellchat_backend(
    *,
    adata,
    tables: Path,
    figures: Path,
    analysis_level: str,
) -> tuple[dict[str, Any], dict[str, str]]:
    backend_id = "scrna.communication.cellchat_optional"
    interactions_path = tables / "cellchat_interactions.tsv"
    pathway_path = tables / "cellchat_pathway_summary.tsv"
    network_path = tables / "cellchat_network_summary.tsv"
    status_path = tables / "cellchat_backend_status.tsv"
    manifest_path = tables / "cellchat_backend_manifest.json"
    script_path = tables / "cellchat_run.R"
    expr_path = tables / "cellchat_input_expression.tsv"
    meta_path = tables / "cellchat_input_metadata.tsv"
    figure_path = figures / "cellchat_network.png"
    artifacts = {
        "cellchat_interactions": str(interactions_path),
        "cellchat_pathway_summary": str(pathway_path),
        "cellchat_network_summary": str(network_path),
        "cellchat_backend_status": str(status_path),
        "cellchat_backend_manifest": str(manifest_path),
        "cellchat_run_script": str(script_path),
        "cellchat_network": str(figure_path),
    }
    rscript = _rscript_str("scrna")
    if not rscript:
        reason = "dependency_missing:Rscript"
        _write_cellchat_skip_outputs(interactions_path, pathway_path, network_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    packages = _available_r_packages(rscript, ("CellChat",))
    if not packages.get("CellChat"):
        reason = "dependency_missing:CellChat"
        _write_cellchat_skip_outputs(interactions_path, pathway_path, network_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    groupby, reason = _communication_groupby(adata)
    if reason:
        _write_cellchat_skip_outputs(interactions_path, pathway_path, network_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    try:
        _write_cellchat_inputs(adata, groupby, expr_path, meta_path)
        _write_cellchat_r_script(script_path, expr_path, meta_path, interactions_path, pathway_path, network_path, status_path, analysis_level)
        completed = subprocess.run([rscript, str(script_path)], check=False, text=True, capture_output=True, timeout=900)
        (tables / "cellchat_stdout.log").write_text(completed.stdout, encoding="utf-8")
        (tables / "cellchat_stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            reason = f"backend_failed:CellChat_returncode_{completed.returncode}"
            _write_cellchat_skip_outputs(interactions_path, pathway_path, network_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
            return _backend_row(backend_id, "failed", analysis_level, reason), artifacts
        if not interactions_path.exists() or interactions_path.stat().st_size == 0:
            reason = "empty_result:cellchat_interactions"
            _write_cellchat_skip_outputs(interactions_path, pathway_path, network_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
            return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
        interactions = pd.read_csv(interactions_path, sep="\t")
        if interactions.empty:
            reason = "empty_result:cellchat_interactions"
            _write_cellchat_skip_outputs(interactions_path, pathway_path, network_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
            return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
        _write_cellchat_network_figure(network_path, figure_path)
        _write_cellchat_manifest(manifest_path, backend_id, "ready", analysis_level, "", artifacts, groupby=groupby)
        return _backend_row(backend_id, "ready", analysis_level, ""), artifacts
    except Exception as exc:
        reason = f"backend_failed:{type(exc).__name__}:{exc}"
        _write_cellchat_skip_outputs(interactions_path, pathway_path, network_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "failed", analysis_level, reason), artifacts


def _run_nichenet_backend(
    *,
    adata,
    tables: Path,
    figures: Path,
    analysis_level: str,
    resource_path: Path | None,
) -> tuple[dict[str, Any], dict[str, str]]:
    backend_id = "scrna.communication.nichenet_optional"
    activity_path = tables / "nichenet_ligand_activity.tsv"
    links_path = tables / "nichenet_ligand_target_links.tsv"
    status_path = tables / "nichenet_backend_status.tsv"
    manifest_path = tables / "nichenet_backend_manifest.json"
    figure_path = figures / "nichenet_activity.png"
    artifacts = {
        "nichenet_ligand_activity": str(activity_path),
        "nichenet_ligand_target_links": str(links_path),
        "nichenet_backend_status": str(status_path),
        "nichenet_backend_manifest": str(manifest_path),
        "nichenet_activity": str(figure_path),
    }
    groupby, reason = _communication_groupby(adata)
    if reason:
        _write_nichenet_skip_outputs(activity_path, links_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    if resource_path is None:
        reason = "missing_ligand_target_resource:provide nichenet_resource or resources.nichenet_ligand_target_matrix"
        _write_nichenet_skip_outputs(activity_path, links_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    if not resource_path.exists():
        reason = f"missing_ligand_target_resource:{resource_path}"
        _write_nichenet_skip_outputs(activity_path, links_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
    try:
        resource = _read_nichenet_resource(resource_path)
        if resource.empty:
            reason = f"empty_ligand_target_resource:{resource_path}"
            _write_nichenet_skip_outputs(activity_path, links_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
            return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
        targets = _nichenet_target_genes(tables / "de_condition.tsv")
        if not targets:
            reason = "missing_target_gene_set:condition_de_or_target_gene_table_required"
            _write_nichenet_skip_outputs(activity_path, links_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
            return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
        activity, links = _score_nichenet_resource(resource, targets, backend_id, analysis_level, groupby, str(resource_path))
        if activity.empty or links.empty:
            reason = "empty_result:no_ligand_target_overlap"
            _write_nichenet_skip_outputs(activity_path, links_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
            return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts
        activity.to_csv(activity_path, sep="\t", index=False)
        links.to_csv(links_path, sep="\t", index=False)
        _write_nichenet_status(status_path, backend_id, analysis_level, "ready", "", resource_path, int(activity.shape[0]), int(links.shape[0]))
        _write_nichenet_manifest(manifest_path, backend_id, "ready", analysis_level, "", artifacts, resource_path=resource_path, groupby=groupby)
        _write_nichenet_activity_figure(activity, figure_path)
        return _backend_row(backend_id, "ready", analysis_level, ""), artifacts
    except Exception as exc:
        reason = f"backend_failed:{type(exc).__name__}:{exc}"
        _write_nichenet_skip_outputs(activity_path, links_path, status_path, manifest_path, figure_path, backend_id, analysis_level, reason, artifacts)
        return _backend_row(backend_id, "failed", analysis_level, reason), artifacts


SCRN_ENRICHMENT_BACKEND_ID = "scrna.enrichment.local_gmt"


def _run_scrna_enrichment(
    *,
    config: dict | None,
    tables: Path,
    output_dir: Path,
    analysis_level: str,
) -> tuple[dict[str, Any] | None, dict[str, str]]:
    """Local-GMT enrichment over scRNA DE/marker tables; explicit skip otherwise."""
    candidate_tables = [
        tables / "pseudobulk_de_results.tsv",
        tables / "de_condition.tsv",
        tables / "marker_genes.tsv",
    ]
    de_table_path = next((path for path in candidate_tables if path.is_file()), None)
    if de_table_path is None:
        return None, {}
    from ultimate.enrichment_backend import run_enrichment

    try:
        de_table = pd.read_csv(de_table_path, sep="\t")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return _backend_row(SCRN_ENRICHMENT_BACKEND_ID, "skipped", analysis_level, "unreadable_de_table"), {}
    result = run_enrichment(
        de_table=de_table,
        config=config if isinstance(config, dict) else {},
        module_name="scrna",
        output_dir=output_dir,
        analysis_level=analysis_level,
    )
    return (
        _backend_row(
            SCRN_ENRICHMENT_BACKEND_ID,
            str(result.get("status") or "not_recorded"),
            analysis_level,
            str(result.get("reason") or ""),
        ),
        {key: str(value) for key, value in (result.get("artifacts") or {}).items()},
    )


PAIRED_PSEUDOBULK_BACKEND_ID = "scrna.pseudobulk.edger_paired"


def _paired_pseudobulk_config(config: dict | None) -> dict[str, Any] | None:
    """Return paired-design keys when the operator enabled paired edgeR DE."""
    cfg = config if isinstance(config, dict) else {}
    modules = cfg.get("modules") if isinstance(cfg.get("modules"), dict) else {}
    scrna_cfg = modules.get("scrna") if isinstance(modules.get("scrna"), dict) else {}
    pseudobulk = scrna_cfg.get("pseudobulk") if isinstance(scrna_cfg.get("pseudobulk"), dict) else {}
    enabled = (
        bool(pseudobulk.get("edger_paired"))
        or bool(pseudobulk.get("paired"))
        or bool(scrna_cfg.get("pseudobulk_paired"))
    )
    donor = pseudobulk.get("donor_column") or scrna_cfg.get("donor_column")
    condition = pseudobulk.get("condition_column") or scrna_cfg.get("condition_column")
    cell_type_col = pseudobulk.get("cell_type_column") or scrna_cfg.get("cell_type_column")
    has_design = bool(donor and condition and cell_type_col)
    if not (enabled or has_design or pseudobulk.get("sce_rda")):
        return None
    return {
        "enabled": enabled or has_design,
        "donor_column": str(donor or "donor"),
        "condition_column": str(condition or "condition"),
        "cell_type_column": str(cell_type_col or "cluster"),
        "cell_type": str(pseudobulk.get("cell_type") or scrna_cfg.get("cell_type") or ""),
        "control": str(pseudobulk.get("control") or scrna_cfg.get("control") or cfg.get("design", {}).get("control") or "control"),
        "case": str(pseudobulk.get("case") or scrna_cfg.get("case") or cfg.get("design", {}).get("case") or "treated"),
        "sce_rda": pseudobulk.get("sce_rda") or scrna_cfg.get("sce_rda"),
        "counts": pseudobulk.get("counts"),
        "design": pseudobulk.get("design"),
        "min_cells": int(pseudobulk.get("min_cells", 1)),
        "min_pairs": int(pseudobulk.get("min_pairs", 4)),
    }


def _run_paired_pseudobulk_edger(
    *,
    config: dict | None,
    tables: Path,
    figures: Path,
    objects: Path,
    analysis_level: str,
) -> tuple[dict[str, Any] | None, dict[str, str]]:
    """Paired edgeR pseudobulk DE inside `ultimate run`, or an explicit skip row."""
    opts = _paired_pseudobulk_config(config)
    if opts is None:
        return None, {}
    backend_id = PAIRED_PSEUDOBULK_BACKEND_ID
    tables.mkdir(parents=True, exist_ok=True)
    status_path = tables / "pseudobulk_paired_backend_status.tsv"
    result_path = tables / "pseudobulk_de_results.tsv"
    artifacts = {
        "pseudobulk_paired_backend_status": str(status_path),
        "pseudobulk_paired_de_results": str(result_path),
    }
    # Route to the R environment that actually owns edgeR.
    if _rscript_str("scrna") == "":
        return _paired_skip(backend_id, analysis_level, "dependency_missing:Rscript", status_path, artifacts)
    if not opts["cell_type"]:
        return _paired_skip(backend_id, analysis_level, "missing_cell_type_for_paired_design", status_path, artifacts)
    if not (opts["sce_rda"] or (opts["counts"] and opts["design"])):
        return _paired_skip(backend_id, analysis_level, "missing_paired_input:sce_rda_or_counts+design", status_path, artifacts)

    from ultimate.scrna_pseudobulk_de import run_scrna_pseudobulk_de

    try:
        manifest = run_scrna_pseudobulk_de(
            output_dir=tables.parent.parent,
            sce_rda=Path(str(opts["sce_rda"])) if opts["sce_rda"] else None,
            counts=Path(str(opts["counts"])) if opts["counts"] else None,
            design=Path(str(opts["design"])) if opts["design"] else None,
            donor_column=opts["donor_column"],
            condition_column=opts["condition_column"],
            cell_type_column=opts["cell_type_column"],
            cell_type=opts["cell_type"],
            control=opts["control"],
            case=opts["case"],
            min_cells=opts["min_cells"],
            min_pairs=opts["min_pairs"],
            analysis_level=analysis_level,
        )
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        return _paired_skip(backend_id, analysis_level, f"paired_pseudobulk_failed:{type(exc).__name__}:{exc}", status_path, artifacts)

    status = str(manifest.get("status") or "not_recorded")
    reason = str(manifest.get("skip_reason") or manifest.get("reason") or "")
    pd.DataFrame(
        [{"backend_id": backend_id, "status": status, "reason": reason, "analysis_level": analysis_level,
          "n_pairs": manifest.get("n_pairs", ""), "engine": manifest.get("engine", "edgeR"),
          "delivery_allowed": False}]
    ).to_csv(status_path, sep="\t", index=False)
    artifacts["pseudobulk_paired_manifest"] = str(tables / "pseudobulk_de_backend_manifest.json")
    row = _backend_row(backend_id, status if status in {"ready", "skipped", "failed"} else "not_recorded", analysis_level, reason)
    return row, artifacts


def _paired_skip(
    backend_id: str,
    analysis_level: str,
    reason: str,
    status_path: Path,
    artifacts: dict[str, str],
) -> tuple[dict[str, Any], dict[str, str]]:
    pd.DataFrame(
        [{"backend_id": backend_id, "status": "skipped", "reason": reason, "analysis_level": analysis_level, "delivery_allowed": False}]
    ).to_csv(status_path, sep="\t", index=False)
    return _backend_row(backend_id, "skipped", analysis_level, reason), artifacts


def _run_pseudobulk_de_backend(*, tables: Path, analysis_level: str) -> tuple[dict[str, Any], dict[str, str]]:
    backend_id = "scrna.pseudobulk.deseq2_edger"
    status_path = tables / "pseudobulk_de_backend_status.tsv"
    result_path = tables / "pseudobulk_de_results.tsv"
    script_path = tables / "pseudobulk_deseq2_edgeR_handoff.R"
    manifest_path = tables / "pseudobulk_de_backend_manifest.json"
    versions_path = tables / "pseudobulk_de_backend_versions.tsv"
    log_path = tables / "pseudobulk_de_backend.log"
    figure_path = tables.parent / "figures" / "pseudobulk_de_volcano.png"
    object_path = tables.parent.parent / "objects" / "scrna_pseudobulk_de_backend.rds"
    counts_path = tables / "pseudobulk_counts.tsv"
    design_path = tables / "pseudobulk_design.tsv"
    reason = _pseudobulk_design_blocker(design_path)
    if not reason:
        rscript = _rscript_str("scrna")
        if not rscript:
            reason = "dependency_missing:Rscript"
        else:
            packages = _available_r_packages(rscript, ("DESeq2", "edgeR"))
            if not (packages.get("DESeq2") or packages.get("edgeR")):
                reason = "dependency_missing:DESeq2_or_edgeR"
            else:
                script = resolve_packaged_asset("scripts", "R", "scrna_pseudobulk_de_backend.R")
                command = [
                    rscript,
                    str(script),
                    "--counts",
                    str(counts_path),
                    "--design",
                    str(design_path),
                    "--tables-dir",
                    str(tables),
                    "--backend-id",
                    backend_id,
                    "--analysis-level",
                    analysis_level,
                ]
                try:
                    completed = subprocess.run(command, check=True, text=True, capture_output=True, timeout=180)
                    log_path.write_text(
                        "COMMAND\t" + " ".join(command) + "\n\nSTDOUT\n" + completed.stdout + "\n\nSTDERR\n" + completed.stderr,
                        encoding="utf-8",
                    )
                    _write_pseudobulk_r_handoff(script_path, counts_path, design_path)
                    return _backend_row(backend_id, "ready", analysis_level, ""), _pseudobulk_artifacts(
                        status_path,
                        result_path,
                        script_path,
                        manifest_path,
                        versions_path,
                        log_path,
                        figure_path,
                        object_path,
                    )
                except subprocess.CalledProcessError as exc:
                    reason = f"backend_failed:Rscript_exit_{exc.returncode}"
                    log_path.write_text(
                        "COMMAND\t" + " ".join(command) + "\n\nSTDOUT\n" + (exc.stdout or "") + "\n\nSTDERR\n" + (exc.stderr or ""),
                        encoding="utf-8",
                    )
                except subprocess.TimeoutExpired as exc:
                    reason = "backend_failed:Rscript_timeout"
                    log_path.write_text(
                        "COMMAND\t" + " ".join(command) + "\n\nSTDOUT\n" + (exc.stdout or "") + "\n\nSTDERR\n" + (exc.stderr or ""),
                        encoding="utf-8",
                    )
                _write_pseudobulk_status(status_path, backend_id, analysis_level, "failed", reason)
                _write_skip_table(result_path, backend_id, analysis_level, reason)
                _write_pseudobulk_backend_manifest(
                    manifest_path,
                    backend_id=backend_id,
                    status="failed",
                    analysis_level=analysis_level,
                    reason=reason,
                    artifacts=_pseudobulk_artifacts(status_path, result_path, script_path, manifest_path, versions_path, log_path),
                )
                _write_pseudobulk_r_handoff(script_path, counts_path, design_path)
                return _backend_row(backend_id, "failed", analysis_level, reason), _pseudobulk_artifacts(
                    status_path,
                    result_path,
                    script_path,
                    manifest_path,
                    versions_path,
                    log_path,
                    )
    _write_pseudobulk_status(status_path, backend_id, analysis_level, "skipped", reason)
    _write_skip_table(result_path, backend_id, analysis_level, reason)
    _write_pseudobulk_backend_manifest(
        manifest_path,
        backend_id=backend_id,
        status="skipped",
        analysis_level=analysis_level,
        reason=reason,
        artifacts=_pseudobulk_artifacts(status_path, result_path, script_path, manifest_path, versions_path, log_path),
    )
    pd.DataFrame([{"package": "DESeq2_or_edgeR", "version": "", "status": "not_run", "reason": reason}]).to_csv(versions_path, sep="\t", index=False)
    log_path.write_text(f"pseudobulk DE backend skipped: {reason}\n", encoding="utf-8")
    _write_pseudobulk_r_handoff(script_path, counts_path, design_path)
    return _backend_row(backend_id, "skipped", analysis_level, reason), _pseudobulk_artifacts(
        status_path,
        result_path,
        script_path,
        manifest_path,
        versions_path,
        log_path,
    )


def _write_figures(sc, adata, figures: Path) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    apply_clinical_journal_style()
    sc.pl.pca(adata, color="condition", show=False)
    save_figure(figures / "pca_condition.png")
    sc.pl.umap(adata, color=["cluster", "condition"], show=False, wspace=0.35)
    save_figure(figures / "umap_cluster_condition.png")

    de = pd.read_csv(figures.parent / "tables" / "de_condition.tsv", sep="\t")
    condition_de_skip_reason = _condition_de_skip_reason(de)
    if condition_de_skip_reason:
        _write_status_figure(figures / "volcano_condition.png", title="Condition DE skipped", message=condition_de_skip_reason)
    elif {"logfoldchanges", "pvals_adj"}.issubset(de.columns):
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        y = -np.log10(np.maximum(de["pvals_adj"].fillna(1).to_numpy(), 1e-300))
        ax.scatter(de["logfoldchanges"].fillna(0), y, s=14, alpha=0.72, color="#607D9B")
        ax.axvline(1, color="#B75E5E", lw=1)
        ax.axvline(-1, color="#4D7EA8", lw=1)
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 adjusted p")
        fig.tight_layout()
        save_figure(figures / "volcano_condition.png")

    comp = pd.read_csv(figures.parent / "tables" / "cell_type_composition.tsv", sep="\t")
    pivot = comp.pivot_table(index="sample_id", columns="cluster", values="fraction", fill_value=0)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    pivot.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    ax.set_ylabel("Fraction")
    ax.set_xlabel("")
    ax.legend(title="Cluster", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    save_figure(figures / "cell_composition.png")

    scores = figures.parent / "tables" / "signature_scores_by_cluster.tsv"
    if scores.exists():
        frame = pd.read_csv(scores, sep="\t").set_index("cluster").select_dtypes(include=[np.number])
        fig, ax = plt.subplots(figsize=(7, max(3.5, 0.35 * len(frame))))
        sns.heatmap(frame, cmap="vlag", center=0, ax=ax)
        fig.tight_layout()
        save_figure(figures / "signature_heatmap.png")


def _add_analysis_level(frame: pd.DataFrame, analysis_level: str) -> pd.DataFrame:
    frame = frame.copy()
    frame["analysis_level"] = analysis_level
    return frame


def _condition_de_blocker(adata) -> str:
    if "condition" not in adata.obs:
        return "missing_condition_column"
    labels = adata.obs["condition"].astype(str).str.strip()
    labels = labels[~labels.isin({"", "nan", "None", "unknown"})]
    counts = labels.value_counts()
    if int(counts.shape[0]) < 2:
        return f"condition_de_requires_at_least_two_conditions:n_conditions={int(counts.shape[0])}"
    if int((counts >= 2).sum()) < 2:
        return "condition_de_requires_two_conditions_with_at_least_two_cells"
    return ""


def _write_condition_de_skip_table(path: Path, adata, analysis_level: str, reason: str) -> None:
    labels = adata.obs["condition"].astype(str).str.strip() if "condition" in adata.obs else pd.Series([], dtype=str)
    labels = labels[~labels.isin({"", "nan", "None", "unknown"})]
    condition_counts = labels.value_counts().sort_index()
    pd.DataFrame(
        [
            {
                "backend_id": "scrna.condition.rank_genes_groups",
                "status": "skipped",
                "reason": reason,
                "n_conditions": int(condition_counts.shape[0]),
                "condition_counts": json.dumps({str(k): int(v) for k, v in condition_counts.items()}, sort_keys=True),
                "analysis_level": analysis_level,
                "delivery_allowed": False,
                "warning": "Condition-level differential expression requires at least two eligible conditions; this table is a status record, not a DE result.",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def _condition_de_skip_reason(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "condition_de_not_available:empty_table"
    if "status" in frame.columns and frame["status"].astype(str).eq("skipped").any():
        skipped = frame[frame["status"].astype(str).eq("skipped")]
        if "reason" in skipped.columns and not skipped["reason"].empty:
            return str(skipped["reason"].iloc[0])
        return "condition_de_skipped"
    if not {"logfoldchanges", "pvals_adj"}.issubset(frame.columns):
        return "condition_de_not_available:missing_logfoldchanges_or_pvals_adj"
    return ""


def _write_annotation_placeholder(adata, tables: Path, analysis_level: str) -> str:
    warning = "cluster labels are computational clusters only; no biological cell type annotation source was provided"
    rows = (
        adata.obs.groupby("cluster", observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .assign(
            cell_type_label=lambda frame: "cluster_" + frame["cluster"].astype(str),
            annotation_status="placeholder_not_cell_type",
            warning=warning,
            analysis_level=analysis_level,
        )
    )
    rows[["cluster", "cell_type_label", "annotation_status", "warning", "n_cells", "analysis_level"]].to_csv(
        tables / "cell_type_annotation_placeholder.tsv",
        sep="\t",
        index=False,
    )
    return warning


def _write_pseudobulk(counts_adata, obs: pd.DataFrame, tables: Path, analysis_level: str) -> dict[str, Path]:
    from scipy import sparse

    obs = obs.loc[counts_adata.obs_names].copy()
    if "sample_id" not in obs:
        obs["sample_id"] = "sample_1"
    if "condition" not in obs:
        obs["condition"] = "unknown"
    if "cluster" not in obs:
        obs["cluster"] = "0"
    obs["sample_id"] = obs["sample_id"].astype(str)
    obs["condition"] = obs["condition"].astype(str)
    obs["cluster"] = obs["cluster"].astype(str)
    groups = obs.groupby(["sample_id", "condition", "cluster"], observed=False).indices
    columns: dict[str, np.ndarray] = {}
    design_rows = []
    matrix = counts_adata.X
    for (sample_id, condition, cluster), indices in groups.items():
        pseudobulk_id = _safe_id(f"{sample_id}__{condition}__cluster{cluster}")
        group_matrix = matrix[list(indices), :]
        summed = group_matrix.sum(axis=0)
        vector = np.asarray(summed).reshape(-1) if sparse.issparse(group_matrix) else np.asarray(summed).reshape(-1)
        columns[pseudobulk_id] = vector
        design_rows.append(
            {
                "pseudobulk_id": pseudobulk_id,
                "sample_id": sample_id,
                "condition": condition,
                "cluster": cluster,
                "n_cells": int(len(indices)),
                "analysis_unit": "sample_condition_cluster",
                "analysis_level": analysis_level,
            }
        )
    counts = pd.DataFrame(columns, index=pd.Index(counts_adata.var_names.astype(str), name="feature_id")).reset_index()
    design = pd.DataFrame(design_rows)
    feature_meta = pd.DataFrame(
        {
            "feature_id": counts_adata.var_names.astype(str),
            "gene_symbol": counts_adata.var_names.astype(str),
            "analysis_level": analysis_level,
        }
    )
    count_path = tables / "pseudobulk_counts.tsv"
    design_path = tables / "pseudobulk_design.tsv"
    feature_path = tables / "pseudobulk_feature_metadata.tsv"
    counts.to_csv(count_path, sep="\t", index=False)
    design.to_csv(design_path, sep="\t", index=False)
    feature_meta.to_csv(feature_path, sep="\t", index=False)
    return {
        "pseudobulk_counts": count_path,
        "pseudobulk_design": design_path,
        "pseudobulk_feature_metadata": feature_path,
    }


def _pseudobulk_design_blocker(design_path: Path) -> str:
    if not design_path.exists():
        return "missing:pseudobulk_design.tsv"
    design = pd.read_csv(design_path, sep="\t")
    if design.empty:
        return "empty:pseudobulk_design.tsv"
    if not {"pseudobulk_id", "condition", "cluster", "sample_id"}.issubset(design.columns):
        return "invalid_design:missing_pseudobulk_condition_cluster_or_sample_id"
    blockers = []
    valid_clusters = 0
    for cluster, frame in design.groupby("cluster", observed=False):
        replicate_counts = frame.groupby("condition")["sample_id"].nunique()
        if replicate_counts.shape[0] < 2:
            blockers.append(f"cluster_{cluster}:need_two_conditions")
        elif int(replicate_counts.min()) < 2:
            blockers.append(f"cluster_{cluster}:need_at_least_two_samples_per_condition")
        else:
            valid_clusters += 1
    if valid_clusters < 1:
        return ";".join(blockers)
    return ""


def _rscript_str(module_name: str) -> str:
    """Rscript for this module, routed through the shared environment resolver."""
    resolved = resolve_rscript(module_name)
    return str(resolved) if resolved is not None else ""


def _available_r_packages(rscript: str, packages: tuple[str, ...]) -> dict[str, bool]:
    script = "cat(paste(installed.packages()[, 'Package'], collapse='\\n'))"
    try:
        completed = subprocess.run([rscript, "-e", script], check=False, text=True, capture_output=True, timeout=45)
    except Exception:
        return {package: False for package in packages}
    installed = set(completed.stdout.splitlines())
    return {package: package in installed for package in packages}


def _write_pseudobulk_status(path: Path, backend_id: str, analysis_level: str, status: str, reason: str) -> None:
    pd.DataFrame(
        [
            {
                "backend_id": backend_id,
                "status": status,
                "reason": reason,
                "analysis_level": analysis_level,
                "warning": "pseudobulk DESeq2/edgeR must use raw counts and sufficient biological replicates",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def _write_pseudobulk_backend_manifest(
    path: Path,
    *,
    backend_id: str,
    status: str,
    analysis_level: str,
    reason: str,
    artifacts: dict[str, str],
) -> None:
    manifest = {
        "backend_id": backend_id,
        "status": status,
        "analysis_level": analysis_level,
        "skip_reason": reason,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "artifacts": artifacts,
        "interpretation_warning": "pseudobulk DESeq2/edgeR uses sample-level pseudobulk counts and must not be interpreted as cell-level independent DE.",
    }
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


def _full_gene_expression_adata(adata):
    if adata.raw is None:
        return adata
    expression = adata.raw.to_adata()
    expression.obs = adata.obs.loc[expression.obs_names].copy()
    expression.uns.update(getattr(adata, "uns", {}))
    return expression


def _cluster_expression_matrix(adata) -> pd.DataFrame:
    from scipy import sparse

    source = _full_gene_expression_adata(adata)
    matrix = source.X
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    frame = pd.DataFrame(values, index=source.obs_names.astype(str), columns=source.var_names.astype(str))
    clusters = adata.obs.loc[source.obs_names, "cluster"].astype(str)
    return frame.groupby(clusters, observed=False).mean()


def _present_gene_sets(features: pd.Index | list[str], gene_sets: dict[str, list[str]], *, min_genes: int) -> dict[str, list[str]]:
    available = set(pd.Index(features).astype(str))
    result: dict[str, list[str]] = {}
    for name, genes in gene_sets.items():
        present = [gene for gene in genes if gene in available]
        if len(present) >= min_genes:
            result[name] = present
    return result


def _gene_set_network(gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for source, genes in gene_sets.items():
        for target in genes:
            rows.append({"source": source, "target": target, "weight": 1.0})
    return pd.DataFrame(rows)


def _activity_result_to_long(
    result,
    *,
    backend_id: str,
    method: str,
    value_name: str,
    analysis_level: str,
) -> pd.DataFrame:
    scores = result[0] if isinstance(result, tuple) else result
    frame = pd.DataFrame(scores)
    frame.index.name = "cluster"
    long = frame.reset_index().melt(id_vars="cluster", var_name="term", value_name=value_name)
    long["backend_id"] = backend_id
    long["method"] = method
    long["analysis_level"] = analysis_level
    long["warning"] = "Expression-derived activity score; not a direct functional assay or mechanism proof."
    return long


def _write_signature_scores_table(tables: Path, output_path: Path, analysis_level: str, backend_id: str) -> None:
    source_path = tables / "signature_scores_by_cluster.tsv"
    if not source_path.exists():
        _write_skip_table(output_path, backend_id, analysis_level, "missing:signature_scores_by_cluster.tsv")
        return
    scores = pd.read_csv(source_path, sep="\t")
    score_cols = [column for column in scores.columns if column.endswith("_score")]
    if not score_cols:
        _write_skip_table(output_path, backend_id, analysis_level, "missing:signature_score_columns")
        return
    long = scores.melt(id_vars=["cluster"], value_vars=score_cols, var_name="signature", value_name="score")
    long["signature"] = long["signature"].str.replace("_score", "", regex=False)
    long["backend_id"] = backend_id
    long["method"] = "scanpy_score_genes_cluster_mean"
    long["analysis_level"] = analysis_level
    long["warning"] = "Signature score is expression-derived and should not be interpreted as functional flux or clinical advice."
    long.to_csv(output_path, sep="\t", index=False)


def _write_functional_status(
    path: Path,
    backend_id: str,
    analysis_level: str,
    status: str,
    reason: str,
    *,
    decoupler_version: str = "",
    gseapy_version: str = "",
) -> None:
    pd.DataFrame(
        [
            {
                "backend_id": backend_id,
                "status": status,
                "reason": reason,
                "decoupler_version": decoupler_version,
                "gseapy_version": gseapy_version,
                "analysis_level": analysis_level,
                "warning": "Signature/pathway/TF activity scores are expression-derived; they are not direct functional assays or mechanism proof.",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def _communication_groupby(adata) -> tuple[str, str]:
    candidates = ("reviewed_cell_type", "cell_type", "manual_cell_type", "celltypist_label", "majority_voting")
    for column in candidates:
        if column not in adata.obs:
            continue
        labels = adata.obs[column].astype(str)
        labels = labels[~labels.isin({"", "nan", "None", "unassigned", "unknown"})]
        if labels.nunique() < 2:
            continue
        counts = labels.value_counts()
        if int((counts >= 10).sum()) < 2:
            return "", f"insufficient_cells_per_group:{column}:need_two_groups_with_min_10_cells"
        return column, ""
    return "", "annotation_required:reviewed_cell_type_or_celltypist_label"


def _write_liana_skip_outputs(
    interactions_path: Path,
    network_path: Path,
    status_path: Path,
    figure_path: Path,
    backend_id: str,
    analysis_level: str,
    reason: str,
) -> None:
    _write_skip_table(interactions_path, backend_id, analysis_level, reason)
    _write_skip_table(network_path, backend_id, analysis_level, reason)
    _write_liana_status(status_path, backend_id, analysis_level, "skipped", reason)
    _write_status_figure(
        figure_path,
        title="LIANA skipped",
        message=reason,
    )


def _liana_network_summary(interactions: pd.DataFrame, backend_id: str, groupby: str, analysis_level: str) -> pd.DataFrame:
    frame = interactions.copy()
    for column in ("lrscore", "magnitude_rank", "specificity_rank"):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    agg: dict[str, tuple[str, str]] = {
        "n_interactions": ("ligand_complex", "count"),
    }
    if "lrscore" in frame:
        agg["mean_lrscore"] = ("lrscore", "mean")
        agg["max_lrscore"] = ("lrscore", "max")
    if "magnitude_rank" in frame:
        agg["mean_magnitude_rank"] = ("magnitude_rank", "mean")
    if "specificity_rank" in frame:
        agg["mean_specificity_rank"] = ("specificity_rank", "mean")
    network = frame.groupby(["source", "target"], observed=False).agg(**agg).reset_index()
    network["backend_id"] = backend_id
    network["method"] = "liana_rank_aggregate_consensus"
    network["groupby"] = groupby
    network["analysis_level"] = analysis_level
    network["warning"] = "Aggregated ligand-receptor candidates; not direct cell-cell contact or causal evidence."
    sort_cols = [column for column in ("mean_magnitude_rank", "mean_specificity_rank", "mean_lrscore") if column in network]
    if sort_cols:
        ascending = [True if "rank" in column else False for column in sort_cols]
        network = network.sort_values(sort_cols, ascending=ascending)
    return network


def _write_liana_status(
    path: Path,
    backend_id: str,
    analysis_level: str,
    status: str,
    reason: str,
    *,
    groupby: str = "",
    liana_version: str = "",
    n_interactions: int = 0,
    n_cell_groups: int = 0,
) -> None:
    pd.DataFrame(
        [
            {
                "backend_id": backend_id,
                "status": status,
                "reason": reason,
                "groupby": groupby,
                "liana_version": liana_version,
                "n_interactions": n_interactions,
                "n_cell_groups": n_cell_groups,
                "analysis_level": analysis_level,
                "warning": "LIANA ligand-receptor results are expression-derived statistical candidates and require biological review.",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def _write_liana_dotplot(interactions: pd.DataFrame, path: Path) -> None:
    import matplotlib.pyplot as plt

    apply_clinical_journal_style()
    frame = interactions.copy()
    if "magnitude_rank" in frame:
        frame["magnitude_rank"] = pd.to_numeric(frame["magnitude_rank"], errors="coerce")
        frame = frame.sort_values("magnitude_rank", ascending=True)
    elif "lrscore" in frame:
        frame["lrscore"] = pd.to_numeric(frame["lrscore"], errors="coerce")
        frame = frame.sort_values("lrscore", ascending=False)
    top = frame.head(25).copy()
    top["pair"] = top["source"].astype(str) + " -> " + top["target"].astype(str)
    top["interaction"] = top["ligand_complex"].astype(str) + " - " + top["receptor_complex"].astype(str)
    size_values = pd.to_numeric(top.get("lrscore", pd.Series([1.0] * len(top))), errors="coerce").fillna(0.5).to_numpy()
    color_values = pd.to_numeric(top.get("magnitude_rank", pd.Series(np.arange(len(top)))), errors="coerce").fillna(0.0).to_numpy()
    fig, ax = plt.subplots(figsize=(8.5, max(4.2, 0.26 * max(1, len(top)))))
    ax.scatter(size_values, np.arange(len(top)), s=np.clip(size_values, 0.05, None) * 90, c=color_values, cmap="viridis_r", alpha=0.82)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["interaction"], fontsize=7)
    ax.set_xlabel("LIANA LR score")
    ax.set_ylabel("")
    ax.set_title("Top ligand-receptor candidates")
    ax.invert_yaxis()
    fig.tight_layout()
    save_figure(path)


def _write_cellchat_inputs(adata, groupby: str, expr_path: Path, meta_path: Path) -> None:
    from scipy import sparse

    source = adata.raw.to_adata() if adata.raw is not None else adata
    matrix = source.X
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    values = np.asarray(values, dtype=float)
    cell_names = source.obs_names.astype(str).to_numpy()
    gene_names = source.var_names.astype(str).to_numpy()
    groups = adata.obs.loc[source.obs_names, groupby].astype(str).to_numpy()
    valid = pd.Series(groups).isin(pd.Series(groups).value_counts()[lambda s: s >= 10].index).to_numpy()
    values = values[valid]
    cell_names = cell_names[valid]
    groups = groups[valid]
    if values.shape[0] > 2000:
        rng = np.random.default_rng(49)
        selected = np.sort(rng.choice(values.shape[0], size=2000, replace=False))
        values = values[selected]
        cell_names = cell_names[selected]
        groups = groups[selected]
    gene_score = values.mean(axis=0)
    selected_genes = np.argsort(gene_score)[::-1][: min(2500, values.shape[1])]
    expr = pd.DataFrame(values[:, selected_genes].T, index=gene_names[selected_genes], columns=cell_names)
    meta = pd.DataFrame({"cell": cell_names, "cell_group": groups})
    expr.index.name = "gene"
    expr.to_csv(expr_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)


def _write_cellchat_r_script(
    script_path: Path,
    expr_path: Path,
    meta_path: Path,
    interactions_path: Path,
    pathway_path: Path,
    network_path: Path,
    status_path: Path,
    analysis_level: str,
) -> None:
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env Rscript",
                "suppressPackageStartupMessages(library(CellChat))",
                f"expr_path <- {json.dumps(str(expr_path))}",
                f"meta_path <- {json.dumps(str(meta_path))}",
                f"interactions_path <- {json.dumps(str(interactions_path))}",
                f"pathway_path <- {json.dumps(str(pathway_path))}",
                f"network_path <- {json.dumps(str(network_path))}",
                f"status_path <- {json.dumps(str(status_path))}",
                f"analysis_level <- {json.dumps(analysis_level)}",
                "expr <- read.delim(expr_path, row.names = 1, check.names = FALSE)",
                "meta <- read.delim(meta_path, check.names = FALSE)",
                "rownames(meta) <- meta$cell",
                "cellchat <- createCellChat(object = as.matrix(expr), meta = meta, group.by = 'cell_group')",
                "db <- tryCatch(CellChatDB.human, error = function(e) NULL)",
                "if (is.null(db)) db <- tryCatch(CellChatDB.mouse, error = function(e) NULL)",
                "if (is.null(db)) stop('CellChatDB.human_or_mouse_missing')",
                "cellchat@DB <- db",
                "cellchat <- subsetData(cellchat)",
                "cellchat <- identifyOverExpressedGenes(cellchat)",
                "cellchat <- identifyOverExpressedInteractions(cellchat)",
                "cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)",
                "cellchat <- filterCommunication(cellchat, min.cells = 10)",
                "cellchat <- computeCommunProbPathway(cellchat)",
                "df.net <- subsetCommunication(cellchat)",
                "if (is.null(df.net) || nrow(df.net) == 0) {",
                "  df.net <- data.frame(source=character(), target=character(), ligand=character(), receptor=character(), prob=numeric())",
                "}",
                "df.net$backend_id <- 'scrna.communication.cellchat_optional'",
                "df.net$analysis_level <- analysis_level",
                "df.net$warning <- 'CellChat ligand-receptor results are expression-derived candidates, not direct mechanism proof.'",
                "write.table(df.net, interactions_path, sep = '\\t', quote = FALSE, row.names = FALSE)",
                "pathways <- tryCatch(cellchat@netP$pathways, error = function(e) character())",
                "pathway <- data.frame(pathway = pathways, backend_id = 'scrna.communication.cellchat_optional', analysis_level = analysis_level, warning = 'Pathway activity is CellChat-inferred and requires review.')",
                "write.table(pathway, pathway_path, sep = '\\t', quote = FALSE, row.names = FALSE)",
                "if (nrow(df.net) > 0 && all(c('source', 'target') %in% colnames(df.net))) {",
                "  network <- aggregate(rep(1, nrow(df.net)), list(source=df.net$source, target=df.net$target), sum)",
                "  colnames(network)[3] <- 'n_interactions'",
                "} else {",
                "  network <- data.frame(source=character(), target=character(), n_interactions=integer())",
                "}",
                "network$backend_id <- 'scrna.communication.cellchat_optional'",
                "network$analysis_level <- analysis_level",
                "network$warning <- 'Aggregated CellChat communication candidates; not causal proof.'",
                "write.table(network, network_path, sep = '\\t', quote = FALSE, row.names = FALSE)",
                "status <- data.frame(backend_id='scrna.communication.cellchat_optional', status='ready', reason='', analysis_level=analysis_level, n_interactions=nrow(df.net), n_groups=length(unique(meta$cell_group)), warning='CellChat requires reviewed labels and biological interpretation.')",
                "write.table(status, status_path, sep = '\\t', quote = FALSE, row.names = FALSE)",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_cellchat_skip_outputs(
    interactions_path: Path,
    pathway_path: Path,
    network_path: Path,
    status_path: Path,
    manifest_path: Path,
    figure_path: Path,
    backend_id: str,
    analysis_level: str,
    reason: str,
    artifacts: dict[str, str],
) -> None:
    for path in (interactions_path, pathway_path, network_path):
        _write_skip_table(path, backend_id, analysis_level, reason)
    _write_cellchat_status(status_path, backend_id, analysis_level, "skipped", reason)
    _write_cellchat_manifest(manifest_path, backend_id, "skipped", analysis_level, reason, artifacts)
    _write_status_figure(figure_path, title="CellChat skipped", message=reason)


def _write_cellchat_status(path: Path, backend_id: str, analysis_level: str, status: str, reason: str) -> None:
    pd.DataFrame(
        [
            {
                "backend_id": backend_id,
                "status": status,
                "reason": reason,
                "analysis_level": analysis_level,
                "warning": "CellChat ligand-receptor results require reviewed cell labels and are not direct mechanism proof.",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def _write_cellchat_manifest(
    path: Path,
    backend_id: str,
    status: str,
    analysis_level: str,
    reason: str,
    artifacts: dict[str, str],
    *,
    groupby: str = "",
) -> None:
    path.write_text(
        json.dumps(
            {
                "backend_id": backend_id,
                "status": status,
                "analysis_level": analysis_level,
                "delivery_allowed": False,
                "validation_evidence_allowed": analysis_level == "validated_backend" and status == "ready",
                "skip_reason": reason,
                "groupby": groupby,
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "artifacts": artifacts,
                "interpretation_warning": "CellChat communication is an expression-derived ligand-receptor inference and must not be reported as mechanism proof.",
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )


def _write_cellchat_network_figure(network_path: Path, figure_path: Path) -> None:
    import matplotlib.pyplot as plt

    apply_clinical_journal_style()
    try:
        network = pd.read_csv(network_path, sep="\t")
    except Exception:
        network = pd.DataFrame()
    if network.empty or "n_interactions" not in network:
        _write_status_figure(figure_path, title="CellChat network", message="No interactions retained")
        return
    top = network.sort_values("n_interactions", ascending=False).head(20).copy()
    labels = top["source"].astype(str) + " -> " + top["target"].astype(str)
    fig, ax = plt.subplots(figsize=(7.5, max(3.0, 0.28 * len(top))))
    ax.barh(labels[::-1], top["n_interactions"].astype(float)[::-1], color="#4777A5")
    ax.set_xlabel("Candidate interactions")
    ax.set_title("CellChat communication network summary")
    fig.tight_layout()
    save_figure(figure_path)


def _write_nichenet_skip_outputs(
    activity_path: Path,
    links_path: Path,
    status_path: Path,
    manifest_path: Path,
    figure_path: Path,
    backend_id: str,
    analysis_level: str,
    reason: str,
    artifacts: dict[str, str],
) -> None:
    _write_nichenet_empty_table(activity_path, backend_id, analysis_level, reason, table_type="ligand_activity")
    _write_nichenet_empty_table(links_path, backend_id, analysis_level, reason, table_type="ligand_target_links")
    _write_nichenet_status(status_path, backend_id, analysis_level, "skipped", reason, None, 0, 0)
    _write_nichenet_manifest(manifest_path, backend_id, "skipped", analysis_level, reason, artifacts, resource_path=None, groupby="")
    _write_status_figure(figure_path, title="NicheNet skipped", message=reason)


def _write_nichenet_empty_table(path: Path, backend_id: str, analysis_level: str, reason: str, *, table_type: str) -> None:
    common = {
        "backend_id": backend_id,
        "status": "skipped",
        "reason": reason,
        "analysis_level": analysis_level,
        "warning": "NicheNet-style ligand-target scoring is a mechanism-hypothesis screen, not mechanism proof.",
    }
    if table_type == "ligand_activity":
        rows = [{"ligand": "", "activity_score": np.nan, "n_target_overlap": 0, **common}]
    else:
        rows = [{"ligand": "", "target": "", "weight": np.nan, "target_in_de": False, **common}]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _read_nichenet_resource(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(path, sep=sep)
    lower = {column.lower(): column for column in frame.columns}
    ligand_col = next((lower[name] for name in ("ligand", "source", "ligand_symbol", "from") if name in lower), None)
    target_col = next((lower[name] for name in ("target", "target_gene", "gene", "target_symbol", "to") if name in lower), None)
    weight_col = next((lower[name] for name in ("weight", "score", "regulatory_potential", "prior") if name in lower), None)
    if ligand_col is None or target_col is None:
        raise ValueError("nichenet resource must contain ligand/source and target/gene columns")
    out = pd.DataFrame(
        {
            "ligand": frame[ligand_col].astype(str).str.strip(),
            "target": frame[target_col].astype(str).str.strip(),
            "weight": pd.to_numeric(frame[weight_col], errors="coerce") if weight_col else 1.0,
        }
    )
    out["weight"] = pd.to_numeric(out["weight"], errors="coerce").fillna(1.0)
    out = out[(out["ligand"] != "") & (out["target"] != "")]
    return out.drop_duplicates(["ligand", "target"]).reset_index(drop=True)


def _nichenet_target_genes(de_path: Path, *, max_targets: int = 250) -> set[str]:
    if not de_path.exists() or de_path.stat().st_size == 0:
        return set()
    frame = pd.read_csv(de_path, sep="\t")
    gene_col = next((column for column in ("names", "gene", "gene_symbol", "feature_id") if column in frame.columns), None)
    if gene_col is None:
        return set()
    ranked = frame.copy()
    if "pvals_adj" in ranked.columns:
        ranked["_rank_p"] = pd.to_numeric(ranked["pvals_adj"], errors="coerce").fillna(1.0)
        ranked = ranked.sort_values("_rank_p", ascending=True)
    elif "pvals" in ranked.columns:
        ranked["_rank_p"] = pd.to_numeric(ranked["pvals"], errors="coerce").fillna(1.0)
        ranked = ranked.sort_values("_rank_p", ascending=True)
    elif "logfoldchanges" in ranked.columns:
        ranked["_rank_lfc"] = pd.to_numeric(ranked["logfoldchanges"], errors="coerce").abs().fillna(0.0)
        ranked = ranked.sort_values("_rank_lfc", ascending=False)
    genes = ranked[gene_col].astype(str).str.strip()
    genes = genes[(genes != "") & (~genes.isin({"nan", "None"}))]
    return set(genes.head(max_targets).str.upper())


def _score_nichenet_resource(
    resource: pd.DataFrame,
    targets: set[str],
    backend_id: str,
    analysis_level: str,
    groupby: str,
    resource_path: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = resource.copy()
    frame["target_upper"] = frame["target"].astype(str).str.upper()
    frame["target_in_de"] = frame["target_upper"].isin(targets)
    links = frame[frame["target_in_de"]].copy()
    if links.empty:
        return pd.DataFrame(), pd.DataFrame()
    links["weight"] = pd.to_numeric(links["weight"], errors="coerce").fillna(1.0)
    activity = (
        links.groupby("ligand", observed=False)
        .agg(activity_score=("weight", "sum"), n_target_overlap=("target", "nunique"), max_weight=("weight", "max"))
        .reset_index()
        .sort_values(["activity_score", "n_target_overlap"], ascending=False)
    )
    activity["rank"] = np.arange(1, len(activity) + 1)
    for table in (activity, links):
        table["backend_id"] = backend_id
        table["status"] = "ready"
        table["analysis_level"] = analysis_level
        table["groupby"] = groupby
        table["resource_path"] = resource_path
        table["method"] = "nichenet_style_ligand_target_overlap_score"
        table["warning"] = "Ligand-target activity is hypothesis generation from expression and prior resources; it is not mechanism proof."
    keep = [
        "ligand",
        "target",
        "weight",
        "target_in_de",
        "backend_id",
        "status",
        "analysis_level",
        "groupby",
        "resource_path",
        "method",
        "warning",
    ]
    return activity, links.sort_values(["ligand", "weight"], ascending=[True, False])[keep].head(500)


def _write_nichenet_status(
    path: Path,
    backend_id: str,
    analysis_level: str,
    status: str,
    reason: str,
    resource_path: Path | None,
    n_ligands: int,
    n_links: int,
) -> None:
    pd.DataFrame(
        [
            {
                "backend_id": backend_id,
                "status": status,
                "reason": reason,
                "resource_path": str(resource_path or ""),
                "n_ligands": n_ligands,
                "n_ligand_target_links": n_links,
                "analysis_level": analysis_level,
                "warning": "NicheNet-style ligand-target analysis is a mechanism-hypothesis screen and requires reviewed sender/receiver labels and target genes.",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def _write_nichenet_manifest(
    path: Path,
    backend_id: str,
    status: str,
    analysis_level: str,
    reason: str,
    artifacts: dict[str, str],
    *,
    resource_path: Path | None,
    groupby: str,
) -> None:
    path.write_text(
        json.dumps(
            {
                "backend_id": backend_id,
                "status": status,
                "analysis_level": analysis_level,
                "delivery_allowed": False,
                "validation_evidence_allowed": analysis_level == "validated_backend" and status == "ready",
                "skip_reason": reason,
                "resource_path": str(resource_path or ""),
                "groupby": groupby,
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "artifacts": artifacts,
                "interpretation_warning": "NicheNet-style ligand-target results generate hypotheses and must not be reported as validated mechanisms.",
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )


def _write_nichenet_activity_figure(activity: pd.DataFrame, figure_path: Path) -> None:
    import matplotlib.pyplot as plt

    apply_clinical_journal_style()
    if activity.empty:
        _write_status_figure(figure_path, title="NicheNet ligand activity", message="No ligand-target overlap")
        return
    top = activity.sort_values("activity_score", ascending=False).head(20).copy()
    fig, ax = plt.subplots(figsize=(7.0, max(3.0, 0.28 * len(top))))
    ax.barh(top["ligand"].astype(str)[::-1], top["activity_score"].astype(float)[::-1], color="#537FA6")
    ax.set_xlabel("Ligand-target overlap score")
    ax.set_title("NicheNet-style ligand activity candidates")
    fig.tight_layout()
    save_figure(figure_path)


def _write_status_figure(path: Path, *, title: str, message: str) -> None:
    import matplotlib.pyplot as plt

    apply_clinical_journal_style()
    fig, ax = plt.subplots(figsize=(6.0, 2.8))
    ax.axis("off")
    ax.text(0.02, 0.72, title, fontsize=13, fontweight="bold", ha="left", va="center")
    ax.text(0.02, 0.42, message[:180], fontsize=9, ha="left", va="center", wrap=True)
    fig.tight_layout()
    save_figure(path)


def _write_pseudobulk_r_handoff(path: Path, counts_path: Path, design_path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env Rscript",
                "# Design-ready handoff for DESeq2/edgeR pseudobulk backend.",
                "# Generated by ultimate validate-scrna; edit contrasts only after reviewing design.",
                f"counts_path <- {json.dumps(str(counts_path))}",
                f"design_path <- {json.dumps(str(design_path))}",
                "counts <- read.delim(counts_path, check.names = FALSE)",
                "design <- read.delim(design_path, check.names = FALSE)",
                "stopifnot('feature_id' %in% colnames(counts))",
                "stopifnot(all(c('pseudobulk_id', 'condition', 'cluster') %in% colnames(design)))",
                "# Guarded backend: run per-cluster DESeq2/edgeR only after replicate, contrast, and raw-count gates pass.",
                "# If any gate fails, this file remains a design-ready handoff and is not a formal DE result.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _pseudobulk_artifacts(
    status_path: Path,
    result_path: Path,
    script_path: Path,
    manifest_path: Path | None = None,
    versions_path: Path | None = None,
    log_path: Path | None = None,
    figure_path: Path | None = None,
    object_path: Path | None = None,
) -> dict[str, str]:
    artifacts = {
        "pseudobulk_de_backend_status": str(status_path),
        "pseudobulk_de_results": str(result_path),
        "pseudobulk_deseq2_edgeR_handoff": str(script_path),
    }
    if manifest_path is not None:
        artifacts["pseudobulk_de_backend_manifest"] = str(manifest_path)
    if versions_path is not None:
        artifacts["pseudobulk_de_backend_versions"] = str(versions_path)
    if log_path is not None:
        artifacts["pseudobulk_de_backend_log"] = str(log_path)
    if figure_path is not None:
        artifacts["pseudobulk_de_volcano"] = str(figure_path)
    if object_path is not None:
        artifacts["pseudobulk_de_rds"] = str(object_path)
    return artifacts


def _celltypist_artifacts(
    annotation_path: Path,
    confidence_path: Path,
    warning_path: Path,
    *,
    decision_path: Path | None = None,
    feature_summary_path: Path | None = None,
    missing_model_path: Path | None = None,
    unmatched_input_path: Path | None = None,
    duplicate_gene_path: Path | None = None,
    label_distribution_path: Path | None = None,
    cluster_crosstab_path: Path | None = None,
    marker_review_path: Path | None = None,
    uncertainty_path: Path | None = None,
    annotated_h5ad_path: Path | None = None,
    manifest_path: Path | None = None,
    report_path: Path | None = None,
    umap_label_path: Path | None = None,
    umap_score_path: Path | None = None,
    marker_dotplot_path: Path | None = None,
) -> dict[str, str]:
    artifacts = {
        "celltypist_annotation": str(annotation_path),
        "annotation_confidence": str(confidence_path),
        "annotation_warning": str(warning_path),
    }
    optional = {
        "celltypist_decision_scores": decision_path,
        "celltypist_feature_match_summary": feature_summary_path,
        "celltypist_model_features_missing_from_input": missing_model_path,
        "celltypist_input_genes_not_in_model": unmatched_input_path,
        "celltypist_duplicate_genes": duplicate_gene_path,
        "celltypist_label_distribution": label_distribution_path,
        "celltypist_cluster_crosstab": cluster_crosstab_path,
        "celltypist_marker_review": marker_review_path,
        "celltypist_uncertain_cells": uncertainty_path,
        "celltypist_annotated_h5ad": annotated_h5ad_path,
        "celltypist_annotation_manifest": manifest_path,
        "celltypist_annotation_report": report_path,
        "celltypist_umap_label": umap_label_path,
        "celltypist_umap_max_probability": umap_score_path,
        "celltypist_marker_dotplot": marker_dotplot_path,
    }
    artifacts.update({key: str(path) for key, path in optional.items() if path is not None})
    return artifacts


def _write_skip_table(path: Path, backend_id: str, analysis_level: str, reason: str) -> None:
    pd.DataFrame(
        [
            {
                "backend_id": backend_id,
                "status": "skipped",
                "reason": reason,
                "analysis_level": analysis_level,
                "delivery_allowed": False,
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def _write_annotation_warning_table(path: Path, backend_id: str, analysis_level: str, warning: str) -> None:
    pd.DataFrame(
        [{"backend_id": backend_id, "annotation_status": "requires_review_or_skipped", "warning": warning, "analysis_level": analysis_level}]
    ).to_csv(path, sep="\t", index=False)


def _backend_row(backend_id: str, status: str, analysis_level: str, reason: str) -> dict[str, Any]:
    return {
        "backend_id": backend_id,
        "status": status,
        "backend_slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "backend_slurm_job_name": os.environ.get("SLURM_JOB_NAME", ""),
        "analysis_level": analysis_level,
        "delivery_allowed": False,
        "validation_evidence_allowed": analysis_level == "validated_backend" and status in {"ready", "design_ready_r_backend_available"},
        "reason": reason,
    }


def _write_backend_execution_manifest(
    *,
    tables: Path,
    backend_rows: list[dict[str, Any]],
    backend_artifacts: dict[str, str],
    analysis_level: str,
) -> Path:
    table_path = tables / "backend_execution.tsv"
    pd.DataFrame(backend_rows).to_csv(table_path, sep="\t", index=False)
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "analysis_level": analysis_level,
        "status_summary": _backend_execution_summary(backend_rows),
        "backend_execution_table": str(table_path),
        "backend_artifacts": backend_artifacts,
        "backends": backend_rows,
        "warning": "skipped/failed backend rows are explicit evidence boundaries and are not formal analysis results",
    }
    path = tables / "backend_execution_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def _backend_execution_summary(backend_rows: list[dict[str, Any]]) -> dict[str, int]:
    summary: dict[str, int] = {}
    for row in backend_rows:
        status = str(row.get("status") or "unknown")
        summary[status] = summary.get(status, 0) + 1
    return summary


def _safe_id(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in value)


def _write_raw_qc_manifest(
    raw_qc: Path,
    *,
    input_path: Path,
    input_type: str,
    tables: Path,
    figures: Path,
    analysis_level: str,
) -> Path:
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "status": "ready",
        "analysis_level": analysis_level,
        "input_path": str(input_path),
        "input_type": input_type,
        "qc_table": str(tables / "qc_metrics.tsv"),
        "qc_figure": str(figures / "qc_violin.png"),
    }
    path = raw_qc / "raw_qc_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def _write_gene_retention_summary(
    path: Path,
    *,
    analysis_level: str,
    n_input_cells: int,
    n_input_genes: int,
    n_qc_cells: int,
    n_qc_genes: int,
    n_hvg_genes: int,
) -> Path:
    pd.DataFrame(
        [
            {
                "stage": "input_after_read_and_optional_sampling",
                "n_cells": int(n_input_cells),
                "n_genes": int(n_input_genes),
                "matrix": "input X",
                "selection_rule": "read input; cells may later be downsampled by --max-cells",
                "analysis_level": analysis_level,
            },
            {
                "stage": "qc_retained_full_gene_space",
                "n_cells": int(n_qc_cells),
                "n_genes": int(n_qc_genes),
                "matrix": "objects/scrna_qc_full_genes.h5ad",
                "selection_rule": "filter_cells(min_genes=min(5,max(1,n_vars//10))) then filter_genes(min_cells=2)",
                "analysis_level": analysis_level,
            },
            {
                "stage": "hvg_analysis_object",
                "n_cells": int(n_qc_cells),
                "n_genes": int(n_hvg_genes),
                "matrix": "objects/scrna_mvp.h5ad",
                "selection_rule": "highly_variable_genes(n_top_genes=min(1000,n_qc_genes), flavor=seurat)",
                "analysis_level": analysis_level,
            },
        ]
    ).to_csv(path, sep="\t", index=False)
    return path


def _write_hvg_genes_table(adata, path: Path, analysis_level: str) -> Path:
    columns = ["highly_variable", "means", "dispersions", "dispersions_norm"]
    frame = pd.DataFrame({"gene": adata.var_names.astype(str)})
    for column in columns:
        if column in adata.var:
            frame[column] = adata.var[column].to_numpy()
    frame["analysis_level"] = analysis_level
    frame.to_csv(path, sep="\t", index=False)
    return path


def _write_matrix_semantics_table(
    path: Path,
    *,
    analysis_level: str,
    hvg_h5ad_path: Path,
    full_gene_h5ad_path: Path,
    n_input_genes: int,
    n_qc_genes: int,
    n_hvg_genes: int,
) -> Path:
    rows = [
        {
            "component": "input_matrix",
            "actual_data_source": "10x/h5ad input read into AnnData.X",
            "gene_range": f"input genes before gene filtering: {int(n_input_genes)}",
            "value_semantics": "raw counts for 10x inputs; h5ad semantics depend on supplied object",
            "method_fit": "source matrix only; not overwritten",
        },
        {
            "component": "objects/scrna_qc_full_genes.h5ad::X",
            "actual_data_source": str(full_gene_h5ad_path),
            "gene_range": f"all QC-retained genes: {int(n_qc_genes)}",
            "value_semantics": "normalize_total(target_sum=1e4)+log1p expression",
            "method_fit": "appropriate for gene-set scoring, annotation-style expression inputs, and full-gene review",
        },
        {
            "component": "objects/scrna_qc_full_genes.h5ad::layers[counts]",
            "actual_data_source": str(full_gene_h5ad_path),
            "gene_range": f"all QC-retained genes: {int(n_qc_genes)}",
            "value_semantics": "raw input counts aligned to the full-gene object",
            "method_fit": "appropriate for count-preserving validation and pseudobulk-style count aggregation",
        },
        {
            "component": "objects/scrna_mvp.h5ad::X",
            "actual_data_source": str(hvg_h5ad_path),
            "gene_range": f"HVG subset: {int(n_hvg_genes)}",
            "value_semantics": "scaled log expression after HVG subsetting",
            "method_fit": "appropriate for PCA/neighbors/UMAP/clustering; not raw counts",
        },
        {
            "component": "objects/scrna_mvp.h5ad::raw.X",
            "actual_data_source": str(hvg_h5ad_path),
            "gene_range": f"all QC-retained genes: {int(n_qc_genes)}",
            "value_semantics": "normalize_total(target_sum=1e4)+log1p expression",
            "method_fit": "used explicitly for marker/condition DE and full-gene expression backends",
        },
        {
            "component": "objects/scrna_mvp.h5ad::layers[counts]",
            "actual_data_source": str(hvg_h5ad_path),
            "gene_range": f"HVG subset: {int(n_hvg_genes)}",
            "value_semantics": "raw input counts aligned to the HVG analysis object",
            "method_fit": "not the full-gene count store; use scrna_qc_full_genes.h5ad for all QC genes",
        },
        {
            "component": "results/tables/marker_genes.tsv",
            "actual_data_source": "scanpy.tl.rank_genes_groups(..., use_raw=True)",
            "gene_range": f"all QC-retained genes through raw.X: {int(n_qc_genes)}",
            "value_semantics": "rank test on full-gene normalized/log expression",
            "method_fit": "descriptive cluster markers, not condition-level replicated DE",
        },
        {
            "component": "results/tables/pseudobulk_counts.tsv",
            "actual_data_source": "counts_adata.X captured before normalization/log1p",
            "gene_range": f"all QC-retained genes: {int(n_qc_genes)}",
            "value_semantics": "sample/condition/cluster sums of raw input counts",
            "method_fit": "appropriate handoff for formal pseudobulk DE gates",
        },
        {
            "component": "CellTypist input when enabled",
            "actual_data_source": "full-gene expression object from adata.raw.to_adata()",
            "gene_range": f"all QC-retained genes: {int(n_qc_genes)}",
            "value_semantics": "normalized/log expression, not scaled HVG X",
            "method_fit": "interface audited; backend remains dependency/model-gated unless configured",
        },
        {
            "component": "LIANA/CellChat/NicheNet expression inputs when enabled",
            "actual_data_source": "adata.raw/full-gene normalized-log expression; CellChat exports up to top expressed genes for R handoff",
            "gene_range": f"all QC-retained genes available before backend-specific caps: {int(n_qc_genes)}",
            "value_semantics": "normalized/log expression for communication-style candidate inference",
            "method_fit": "interface audited; optional backends remain dependency/resource-gated unless configured",
        },
    ]
    frame = pd.DataFrame(rows)
    frame["analysis_level"] = analysis_level
    frame.to_csv(path, sep="\t", index=False)
    return path


def _reproducible_command(
    *,
    input_path: Path,
    input_type: str,
    output_dir: Path,
    samplesheet: Path | None,
    max_cells: int,
    random_seed: int,
    analysis_level: str,
    public_dataset: bool,
    dataset_label: str | None,
    production_approval: dict[str, Any] | None,
    celltypist_model: Path | None = None,
    nichenet_resource: Path | None = None,
) -> str:
    pieces = [
        "ultimate",
        "validate-scrna",
        "--input-path",
        str(input_path),
        "--input-type",
        input_type,
        "--output-dir",
        str(output_dir),
        "--max-cells",
        str(max_cells),
        "--random-seed",
        str(random_seed),
        "--analysis-level",
        analysis_level,
    ]
    if samplesheet is not None:
        pieces.extend(["--samplesheet", str(samplesheet)])
    if public_dataset:
        pieces.append("--public-dataset")
    if dataset_label:
        pieces.extend(["--dataset-label", str(dataset_label)])
    approval_path = (production_approval or {}).get("_approval_path")
    if approval_path:
        pieces.extend(["--production-approval", str(approval_path)])
    if celltypist_model is not None:
        pieces.extend(["--celltypist-model", str(celltypist_model)])
    if nichenet_resource is not None:
        pieces.extend(["--nichenet-resource", str(nichenet_resource)])
    return " ".join(pieces)


def _approval_summary(approval: dict[str, Any]) -> dict[str, Any]:
    return {
        "approved": bool(approval.get("approved")),
        "approved_by": str(approval.get("approved_by", "")),
        "approved_at": str(approval.get("approved_at", "")),
        "project_id": str(approval.get("project_id", "")),
        "input_path": str(approval.get("input_path", "")),
        "output_dir": str(approval.get("output_dir", "")),
        "delivery_scope": str(approval.get("delivery_scope", "")),
        "reason": str(approval.get("reason", "")),
        "approval_path": str(approval.get("_approval_path", "")),
    }


def _figure_manifest(figures: Path, output_path: Path, analysis_level: str) -> Path:
    rows = [
        {
            "figure": path.name,
            "path": str(path),
            "format": path.suffix.lstrip("."),
            "style": "ultimate_active_style",
            "analysis_level": analysis_level,
        }
        for path in sorted(figures.glob("*.png"))
    ]
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
    return output_path


def _write_report(manifest: dict[str, Any], md_path: Path, html_path: Path) -> None:
    md = [
        "# scRNA-seq MVP 验证报告",
        "",
        f"输入类型：`{manifest['input_type']}`",
        f"输入路径：`{manifest['input_path']}`",
        f"状态：`{manifest['status']}`",
        f"analysis_level：`{manifest['analysis_level']}`",
        f"delivery_allowed：`{manifest['delivery_allowed']}`",
        f"non_delivery_reason：`{manifest['non_delivery_reason'] or 'none'}`",
        "",
        "## 交付边界",
        f"- 当前结果级别：`{manifest['analysis_level']}`",
        f"- 是否允许作为客户正式交付：`{manifest['delivery_allowed']}`",
        f"- 说明：{manifest['non_delivery_reason'] or 'production_backend 可作为正式交付结果。'}",
        f"- 细胞类型注释警示：{manifest['cell_type_annotation_warning']}",
        f"- backend 执行摘要：`{json.dumps(manifest.get('backend_execution_status', {}), ensure_ascii=False)}`",
        "",
        "## 数据概览",
        f"- 细胞数：{manifest['n_cells']}",
        f"- HVG 分析对象基因数：{manifest['n_genes']}",
        f"- 输入基因数：{manifest.get('n_input_genes', 'not_recorded')}",
        f"- QC 后全基因数：{manifest.get('n_qc_genes', 'not_recorded')}",
        f"- HVG 基因数：{manifest.get('n_hvg_genes', manifest['n_genes'])}",
        f"- 全基因对象：`{manifest.get('objects', {}).get('full_gene_h5ad', 'not_recorded')}`",
        f"- 矩阵语义表：`{manifest.get('matrix_semantics', 'not_recorded')}`",
        "",
        "## 已验证步骤",
        "- 输入读取、QC、过滤、归一化、高变基因、PCA、UMAP、聚类、marker gene、条件差异、细胞组成、基础 signature 富集、pseudobulk design-ready matrix、h5ad 导出、图表和 manifest。",
        "- cell type annotation 当前是 cluster placeholder，不能当作真实细胞类型结论。",
        "- Scrublet、CellTypist、LIANA、pseudobulk DESeq2/edgeR 等 backend 会在具备依赖和输入条件时执行；缺依赖、缺模型、缺人工/参考注释或设计不足会写入 backend_execution_manifest，而不是伪装成完成。",
        "- LIANA 细胞通讯结果是基于表达量的配体-受体候选互作推断，不是直接细胞接触或机制证明。",
        "- pseudobulk DE 当前写出真实 R handoff 和设计检查；只有 Rscript 与 DESeq2/edgeR、重复数和 contrasts 齐全后才可升级为正式 DE 结果。",
        "",
        "## 复现命令",
        f"`{manifest['reproducible_command']}`",
    ]
    md_path.write_text("\n".join(md) + "\n", encoding="utf-8")
    html_path.write_text("<html><body>" + "".join(f"<p>{line}</p>" for line in md) + "</body></html>", encoding="utf-8")
