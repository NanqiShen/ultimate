from __future__ import annotations

import json
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import types

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import ultimate.scrna_smoke as scrna_smoke
from ultimate.cli import main
from ultimate.config import dump_yaml
from ultimate.pipeline import run_pipeline_from_config
from ultimate.scrna_smoke import create_demo_inputs, run_scrna_validation
from ultimate.scrna_velocity_backend import has_scrna_velocity_backend_config


def test_create_scrna_demo_inputs(tmp_path: Path) -> None:
    pytest.importorskip("h5py")
    pytest.importorskip("scipy")
    manifest = create_demo_inputs(tmp_path / "demo", n_cells=16, n_genes=32, seed=3)

    assert Path(manifest["tenx_h5"]).exists()
    assert (Path(manifest["tenx_mtx"]) / "matrix.mtx.gz").exists()
    assert (Path(manifest["tenx_mtx"]) / "barcodes.tsv.gz").exists()
    assert (Path(manifest["tenx_mtx"]) / "features.tsv.gz").exists()
    assert Path(manifest["samplesheet"]).exists()
    assert manifest["n_cells"] == 16
    assert manifest["n_genes"] >= 32


def test_create_scrna_demo_inputs_cli(tmp_path: Path) -> None:
    if importlib.util.find_spec("h5py") is None or importlib.util.find_spec("scipy") is None:
        pytest.skip("h5py and scipy are required to materialize 10x demo files")
    runner = CliRunner()
    result = runner.invoke(main, ["create-scrna-demo-inputs", "--output-dir", str(tmp_path / "demo"), "--n-cells", "10", "--n-genes", "30"])

    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert Path(payload["tenx_h5"]).exists()
    assert (Path(payload["tenx_mtx"]) / "matrix.mtx.gz").exists()


def test_scrublet_backend_skips_inputs_too_small_for_pca(tmp_path: Path) -> None:
    class SmallAdata:
        n_obs = 16
        n_vars = 12
        shape = (16, 12)

    tables = tmp_path / "tables"
    figures = tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()

    row, artifacts = scrna_smoke._run_scrublet_backend(
        counts_adata=SmallAdata(),
        adata=SmallAdata(),
        tables=tables,
        figures=figures,
        analysis_level="production_backend",
        random_seed=1,
    )

    assert row["status"] == "skipped"
    assert row["reason"].startswith("input_too_small_for_scrublet")
    assert "backend_failed" not in row["reason"]
    assert Path(artifacts["doublet_scores"]).exists()
    assert Path(artifacts["doublet_summary"]).exists()


def test_scrublet_backend_retries_with_fewer_principal_components(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "scrublet":
            return object()
        return original_find_spec(name)

    calls: list[int] = []

    class FakeScrublet:
        def __init__(self, matrix, random_state=None):
            self.matrix = matrix
            self.random_state = random_state

        def scrub_doublets(self, *, n_prin_comps: int, min_gene_variability_pctl: float = 85, verbose: bool = False):
            calls.append(n_prin_comps)
            if n_prin_comps > 10:
                raise ValueError("n_components=30 must be between 1 and min(n_samples, n_features)=10 with svd_solver='arpack'")
            return np.linspace(0.0, 1.0, self.matrix.shape[0]), np.zeros(self.matrix.shape[0], dtype=bool)

    fake_module = types.SimpleNamespace(Scrublet=FakeScrublet)
    monkeypatch.setitem(sys.modules, "scrublet", fake_module)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)

    class Adata:
        def __init__(self):
            self.X = np.ones((64, 64), dtype=float)
            self.n_obs = 64
            self.n_vars = 64
            self.shape = self.X.shape
            self.obs = pd.DataFrame(index=[f"cell{i}" for i in range(64)])
            self.obs_names = self.obs.index

    tables = tmp_path / "tables"
    figures = tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()
    adata = Adata()
    row, artifacts = scrna_smoke._run_scrublet_backend(
        counts_adata=adata,
        adata=adata,
        tables=tables,
        figures=figures,
        analysis_level="smoke_backend",
        random_seed=1,
    )

    assert row["status"] == "ready"
    summary = pd.read_csv(artifacts["doublet_summary"], sep="\t")
    assert calls == [30, 20, 15, 10]
    assert int(summary.loc[0, "n_prin_comps"]) == 10
    assert float(summary.loc[0, "min_gene_variability_pctl"]) == 85
    assert adata.obs["scrublet_score"].notna().all()


def test_scrublet_backend_retries_with_lower_gene_variability_percentile(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "scrublet":
            return object()
        return original_find_spec(name)

    calls: list[tuple[float, int]] = []

    class FakeScrublet:
        def __init__(self, matrix, random_state=None):
            self.matrix = matrix
            self.random_state = random_state

        def scrub_doublets(self, *, n_prin_comps: int, min_gene_variability_pctl: float = 85, verbose: bool = False):
            calls.append((min_gene_variability_pctl, n_prin_comps))
            if min_gene_variability_pctl > 25:
                raise ValueError("Found array with 0 feature(s) (shape=(64, 0)) while a minimum of 1 is required by PCA.")
            if n_prin_comps > 5:
                raise ValueError("n_components=10 must be strictly less than min(n_samples, n_features)=10 with svd_solver='arpack'")
            return np.linspace(0.0, 1.0, self.matrix.shape[0]), np.zeros(self.matrix.shape[0], dtype=bool)

    fake_module = types.SimpleNamespace(Scrublet=FakeScrublet)
    monkeypatch.setitem(sys.modules, "scrublet", fake_module)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)

    class Adata:
        def __init__(self):
            self.X = np.ones((64, 64), dtype=float)
            self.n_obs = 64
            self.n_vars = 64
            self.shape = self.X.shape
            self.obs = pd.DataFrame(index=[f"cell{i}" for i in range(64)])
            self.obs_names = self.obs.index

    tables = tmp_path / "tables"
    figures = tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()
    adata = Adata()
    row, artifacts = scrna_smoke._run_scrublet_backend(
        counts_adata=adata,
        adata=adata,
        tables=tables,
        figures=figures,
        analysis_level="smoke_backend",
        random_seed=1,
    )

    assert row["status"] == "ready"
    summary = pd.read_csv(artifacts["doublet_summary"], sep="\t")
    assert calls[-1] == (25, 5)
    assert int(summary.loc[0, "n_prin_comps"]) == 5
    assert float(summary.loc[0, "min_gene_variability_pctl"]) == 25
    assert adata.obs["scrublet_score"].notna().all()


def test_validate_scrna_h5ad_outputs_scrna_mvp_artifacts(tmp_path: Path) -> None:
    _require_scrna_runtime()
    demo = create_demo_inputs(tmp_path / "demo_h5ad", n_cells=48, n_genes=60, seed=11)
    manifest = run_scrna_validation(
        input_path=Path(demo["h5ad"]),
        input_type="h5ad",
        output_dir=tmp_path / "run_h5ad",
        samplesheet=Path(demo["samplesheet"]),
        max_cells=48,
    )
    _assert_mvp_outputs(manifest)
    assert manifest["analysis_level"] == "demo_result"
    assert manifest["delivery_allowed"] is False
    assert manifest["non_delivery_reason"] == "generated_demo_data_not_customer_delivery"


def test_validate_scrna_h5ad_prefers_counts_layer_over_logged_x(tmp_path: Path) -> None:
    _require_scrna_runtime()
    ad = pytest.importorskip("anndata")
    sc = pytest.importorskip("scanpy")
    sparse = pytest.importorskip("scipy.sparse")
    rng = np.random.default_rng(31)
    counts = rng.poisson(lam=2.0, size=(36, 72)).astype(float)
    counts[:2, :] += 1.0
    obs_names = [f"cell{i}" for i in range(counts.shape[0])]
    var_names = [f"Gene{i}" for i in range(counts.shape[1])]
    source = ad.AnnData(
        X=sparse.csr_matrix(counts),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    logged = source.copy()
    sc.pp.normalize_total(logged, target_sum=1e4)
    sc.pp.log1p(logged)
    logged.layers["counts"] = source.X.copy()
    input_path = tmp_path / "logged_x_with_counts_layer.h5ad"
    logged.write_h5ad(input_path)
    samplesheet = tmp_path / "samplesheet.tsv"
    pd.DataFrame(
        {
            "barcode": obs_names,
            "sample_id": ["S1"] * len(obs_names),
            "condition": ["only_condition"] * len(obs_names),
        }
    ).to_csv(samplesheet, sep="\t", index=False)

    manifest = run_scrna_validation(
        input_path=input_path,
        input_type="h5ad",
        output_dir=tmp_path / "run_h5ad_counts_layer",
        samplesheet=samplesheet,
        max_cells=36,
    )

    full = sc.read_h5ad(manifest["objects"]["full_gene_h5ad"])
    expected = source.copy()
    sc.pp.filter_cells(expected, min_genes=min(5, max(1, expected.n_vars // 10)))
    sc.pp.filter_genes(expected, min_cells=2)
    assert manifest["pseudobulk_matrix_source"] == "raw_counts_or_input_counts"
    assert full.uns["ultimate_matrix_semantics"]["input_type"] == "h5ad"
    assert full.uns["ultimate_input_matrix_source"] == "h5ad.layers[counts]"
    count_diff = full.layers["counts"] - expected[full.obs_names, full.var_names].X
    count_max_abs = 0.0 if count_diff.nnz == 0 else float(np.max(np.abs(count_diff.data)))
    assert count_max_abs == 0.0

    expected_log = expected[full.obs_names, full.var_names].copy()
    sc.pp.normalize_total(expected_log, target_sum=1e4)
    sc.pp.log1p(expected_log)
    log_diff = full.X - expected_log.X
    log_max_abs = 0.0 if sparse.issparse(log_diff) and log_diff.nnz == 0 else float(np.max(np.abs(log_diff.data if sparse.issparse(log_diff) else np.asarray(log_diff))))
    assert log_max_abs < 1e-6


def test_validate_scrna_h5ad_rejects_logged_x_without_counts_layer(tmp_path: Path) -> None:
    _require_scrna_runtime()
    ad = pytest.importorskip("anndata")
    sc = pytest.importorskip("scanpy")
    sparse = pytest.importorskip("scipy.sparse")
    counts = sparse.csr_matrix(np.array([[1, 0, 2], [0, 3, 1], [2, 1, 0]], dtype=float))
    logged = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(index=["cell1", "cell2", "cell3"]),
        var=pd.DataFrame(index=["Gene1", "Gene2", "Gene3"]),
    )
    sc.pp.normalize_total(logged, target_sum=1e4)
    sc.pp.log1p(logged)
    input_path = tmp_path / "logged_without_counts_layer.h5ad"
    logged.write_h5ad(input_path)

    with pytest.raises(ValueError, match="provide raw counts"):
        run_scrna_validation(
            input_path=input_path,
            input_type="h5ad",
            output_dir=tmp_path / "run_h5ad_reject",
            max_cells=3,
        )


def test_validate_scrna_10x_h5_outputs_scrna_mvp_artifacts(tmp_path: Path) -> None:
    _require_scrna_runtime()
    demo = create_demo_inputs(tmp_path / "demo_h5", n_cells=48, n_genes=60, seed=12)
    manifest = run_scrna_validation(
        input_path=Path(demo["tenx_h5"]),
        input_type="10x_h5",
        output_dir=tmp_path / "run_10x_h5",
        samplesheet=Path(demo["samplesheet"]),
        max_cells=48,
    )
    _assert_mvp_outputs(manifest)
    assert manifest["analysis_level"] == "demo_result"


def test_validate_scrna_10x_mtx_outputs_scrna_mvp_artifacts(tmp_path: Path) -> None:
    _require_scrna_runtime()
    demo = create_demo_inputs(tmp_path / "demo_mtx", n_cells=48, n_genes=60, seed=13)
    manifest = run_scrna_validation(
        input_path=Path(demo["tenx_mtx"]),
        input_type="10x_mtx",
        output_dir=tmp_path / "run_10x_mtx",
        samplesheet=Path(demo["samplesheet"]),
        max_cells=48,
    )
    _assert_mvp_outputs(manifest)
    assert manifest["analysis_level"] == "demo_result"


def test_validate_scrna_single_condition_records_condition_de_skip(tmp_path: Path) -> None:
    _require_scrna_runtime()
    demo = create_demo_inputs(tmp_path / "demo_single_condition", n_cells=48, n_genes=60, seed=16)
    samples = pd.read_csv(demo["samplesheet"], sep="\t")
    samples["condition"] = "healthy_donor"
    samplesheet = tmp_path / "single_condition_samples.tsv"
    samples.to_csv(samplesheet, sep="\t", index=False)

    manifest = run_scrna_validation(
        input_path=Path(demo["tenx_mtx"]),
        input_type="10x_mtx",
        output_dir=tmp_path / "run_single_condition",
        samplesheet=samplesheet,
        max_cells=48,
    )

    _assert_mvp_outputs(manifest)
    de = pd.read_csv(Path(manifest["output_dir"]) / "results" / "tables" / "de_condition.tsv", sep="\t")
    assert de.loc[0, "status"] == "skipped"
    assert "condition_de_requires_at_least_two_conditions" in de.loc[0, "reason"]
    assert int(de.loc[0, "n_conditions"]) == 1
    assert "logfoldchanges" not in de.columns
    volcano = Path(manifest["output_dir"]) / "results" / "figures" / "volcano_condition.png"
    assert volcano.exists()
    assert volcano.stat().st_size > 0


def test_validate_scrna_preserves_full_qc_gene_space_and_raw_counts(tmp_path: Path) -> None:
    _require_scrna_runtime()
    sc = pytest.importorskip("scanpy")
    sparse = pytest.importorskip("scipy.sparse")
    demo = create_demo_inputs(tmp_path / "demo_full_gene_retention", n_cells=80, n_genes=1200, seed=23)

    manifest = run_scrna_validation(
        input_path=Path(demo["tenx_mtx"]),
        input_type="10x_mtx",
        output_dir=tmp_path / "run_full_gene_retention",
        samplesheet=Path(demo["samplesheet"]),
        max_cells=80,
    )

    run_dir = Path(manifest["output_dir"])
    full = sc.read_h5ad(manifest["objects"]["full_gene_h5ad"])
    hvg = sc.read_h5ad(manifest["objects"]["h5ad"])
    source = sc.read_10x_mtx(demo["tenx_mtx"], var_names="gene_symbols", cache=False)
    source.var_names_make_unique()

    assert full.n_obs == hvg.n_obs
    assert full.n_vars > hvg.n_vars
    assert int(manifest["n_qc_genes"]) == full.n_vars
    assert int(manifest["n_hvg_genes"]) == hvg.n_vars
    assert "counts" in full.layers
    assert full.layers["counts"].shape == full.shape

    non_hvg = full.var_names[~full.var["highly_variable"].astype(bool)][0]
    assert non_hvg not in set(hvg.var_names)
    assert non_hvg in set(full.var_names)
    assert non_hvg in set(hvg.raw.var_names)

    expected = source[full.obs_names.astype(str), full.var_names.astype(str)].X
    diff = full.layers["counts"] - expected
    max_abs = 0.0 if diff.nnz == 0 else float(np.max(np.abs(diff.data)))
    assert max_abs == 0.0

    expected_log = source[full.obs_names.astype(str), full.var_names.astype(str)].copy()
    sc.pp.normalize_total(expected_log, target_sum=1e4)
    sc.pp.log1p(expected_log)
    log_diff = full.X - expected_log.X
    log_max_abs = 0.0 if sparse.issparse(log_diff) and log_diff.nnz == 0 else float(np.max(np.abs(log_diff.data if sparse.issparse(log_diff) else np.asarray(log_diff))))
    assert log_max_abs < 1e-6

    marker_names = set(pd.read_csv(run_dir / "results" / "tables" / "marker_genes.tsv", sep="\t")["names"].astype(str))
    assert non_hvg in marker_names

    pseudobulk = pd.read_csv(run_dir / "results" / "tables" / "pseudobulk_counts.tsv", sep="\t").set_index("feature_id")
    assert non_hvg in pseudobulk.index
    design = pd.read_csv(run_dir / "results" / "tables" / "pseudobulk_design.tsv", sep="\t")
    row = design.iloc[0]
    cells = full.obs.index[
        (full.obs["sample_id"].astype(str) == str(row["sample_id"]))
        & (full.obs["condition"].astype(str) == str(row["condition"]))
        & (full.obs["cluster"].astype(str) == str(row["cluster"]))
    ]
    gene_idx = full.var_names.get_loc(non_hvg)
    manual_sum = full[cells, :].layers["counts"][:, gene_idx].sum()
    assert float(pseudobulk.loc[non_hvg, row["pseudobulk_id"]]) == float(manual_sum)

    retention = pd.read_csv(run_dir / "results" / "tables" / "gene_retention_summary.tsv", sep="\t")
    assert {"input_after_read_and_optional_sampling", "qc_retained_full_gene_space", "hvg_analysis_object"} <= set(retention["stage"])
    semantics = pd.read_csv(run_dir / "results" / "tables" / "matrix_semantics.tsv", sep="\t")
    assert "objects/scrna_qc_full_genes.h5ad::layers[counts]" in set(semantics["component"])


def test_celltypist_backend_uses_full_gene_raw_expression(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ad = pytest.importorskip("anndata")

    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "celltypist":
            return object()
        return original_find_spec(name)

    observed: dict[str, object] = {}

    class FakeCellTypist:
        class models:
            class Model:
                @staticmethod
                def load(path: str):
                    return types.SimpleNamespace(features=["G1", "G2", "G3_non_hvg", "G4_non_hvg"])

        @staticmethod
        def annotate(adata, *, model, majority_voting: bool, over_clustering=None):
            observed["shape"] = adata.shape
            observed["genes"] = list(adata.var_names.astype(str))
            observed["min_x"] = float(np.min(np.asarray(adata.X)))
            observed["over_clustering"] = over_clustering
            labels = pd.DataFrame(
                {"predicted_labels": ["mock"] * adata.n_obs, "majority_voting": ["mock_majority"] * adata.n_obs},
                index=adata.obs_names,
            )
            probabilities = pd.DataFrame({"mock": [0.91] * adata.n_obs, "mock_majority": [0.82] * adata.n_obs}, index=adata.obs_names)
            decisions = pd.DataFrame({"mock": [2.2] * adata.n_obs, "mock_majority": [1.5] * adata.n_obs}, index=adata.obs_names)
            return types.SimpleNamespace(predicted_labels=labels, probability_matrix=probabilities, decision_matrix=decisions)

    full = ad.AnnData(
        X=np.log1p(np.array([[1, 2, 3, 4], [2, 1, 0, 5], [3, 0, 1, 6]], dtype=float)),
        obs=pd.DataFrame(index=["cell1", "cell2", "cell3"]),
        var=pd.DataFrame(index=["G1", "G2", "G3_non_hvg", "G4_non_hvg"]),
    )
    hvg = full[:, ["G1", "G2"]].copy()
    hvg.raw = full
    hvg.X = np.array([[-1.0, 0.0], [0.5, -0.2], [1.0, 0.3]])
    hvg.obs["cluster"] = ["0", "0", "1"]

    monkeypatch.setitem(sys.modules, "celltypist", FakeCellTypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", FakeCellTypist.models)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)
    model = tmp_path / "mock.pkl"
    model.write_text("mock", encoding="utf-8")
    tables = tmp_path / "tables"
    tables.mkdir()

    row, _ = scrna_smoke._run_celltypist_backend(adata=hvg, tables=tables, analysis_level="validated_backend", model_path=model)

    assert row["status"] == "ready"
    assert observed["shape"] == (3, 4)
    assert observed["genes"] == ["G1", "G2", "G3_non_hvg", "G4_non_hvg"]
    assert observed["min_x"] >= 0.0
    assert observed["over_clustering"] == "cluster"
    annotation = pd.read_csv(tables / "celltypist_annotation.tsv", sep="\t")
    assert {"celltypist_raw_label", "celltypist_majority_label", "celltypist_max_probability", "uncertain"}.issubset(annotation.columns)
    assert set(annotation["celltypist_raw_label"]) == {"mock"}
    assert set(annotation["celltypist_majority_label"]) == {"mock_majority"}
    assert annotation["uncertain"].eq(False).all()
    assert (tables / "celltypist_feature_match_summary.tsv").exists()
    assert (tables / "celltypist_decision_scores.tsv").exists()


def test_celltypist_backend_aligns_predictions_by_cell_id_and_retains_uncertainty(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ad = pytest.importorskip("anndata")

    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "celltypist":
            return object()
        return original_find_spec(name)

    class FakeCellTypist:
        class models:
            class Model:
                @staticmethod
                def load(path: str):
                    return types.SimpleNamespace(features=["G1", "G2", "G3"])

        @staticmethod
        def annotate(adata, *, model, majority_voting: bool, over_clustering=None):
            shuffled = ["cell3", "cell1", "cell2"]
            labels = pd.DataFrame(
                {
                    "predicted_labels": ["B cells", "T cells", "Monocytes"],
                    "majority_voting": ["B cells", "T cells", "Monocytes"],
                },
                index=shuffled,
            )
            probabilities = pd.DataFrame(
                {
                    "B cells": [0.91, 0.10, 0.20],
                    "T cells": [0.05, 0.93, 0.40],
                    "Monocytes": [0.04, 0.02, 0.45],
                },
                index=shuffled,
            )
            decisions = pd.DataFrame(
                {
                    "B cells": [2.0, -1.0, -0.5],
                    "T cells": [-2.0, 2.1, -0.2],
                    "Monocytes": [-2.1, -2.0, -0.1],
                },
                index=shuffled,
            )
            return types.SimpleNamespace(predicted_labels=labels, probability_matrix=probabilities, decision_matrix=decisions)

    monkeypatch.setitem(sys.modules, "celltypist", FakeCellTypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", FakeCellTypist.models)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)
    adata = ad.AnnData(
        X=np.log1p(np.array([[1, 2, 3], [2, 0, 1], [0, 3, 2]], dtype=float)),
        obs=pd.DataFrame({"cluster": ["0", "1", "0"]}, index=["cell1", "cell2", "cell3"]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    model = tmp_path / "mock.pkl"
    model.write_text("mock", encoding="utf-8")
    tables = tmp_path / "tables"
    objects = tmp_path / "objects"
    manifests = tmp_path / "manifests"
    for directory in (tables, objects, manifests):
        directory.mkdir()

    row, artifacts = scrna_smoke._run_celltypist_backend(
        adata=adata,
        tables=tables,
        objects=objects,
        manifests=manifests,
        analysis_level="validated_backend",
        model_path=model,
        uncertainty_threshold=0.5,
    )

    assert row["status"] == "ready"
    annotation = pd.read_csv(artifacts["celltypist_annotation"], sep="\t")
    assert annotation["barcode"].tolist() == ["cell1", "cell2", "cell3"]
    assert annotation["celltypist_raw_label"].tolist() == ["T cells", "Monocytes", "B cells"]
    assert annotation["celltypist_max_probability"].round(2).tolist() == [0.93, 0.45, 0.91]
    assert annotation["uncertain"].tolist() == [False, True, False]
    uncertain = pd.read_csv(artifacts["celltypist_uncertain_cells"], sep="\t")
    assert uncertain["barcode"].tolist() == ["cell2"]
    assert Path(artifacts["celltypist_annotated_h5ad"]).exists()


def test_celltypist_backend_fails_when_no_model_genes_match(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ad = pytest.importorskip("anndata")

    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "celltypist":
            return object()
        return original_find_spec(name)

    class FakeCellTypist:
        class models:
            class Model:
                @staticmethod
                def load(path: str):
                    return types.SimpleNamespace(features=["MissingA", "MissingB"])

        @staticmethod
        def annotate(adata, *, model, majority_voting: bool, over_clustering=None):
            raise AssertionError("annotate should not be called when feature intersection is empty")

    monkeypatch.setitem(sys.modules, "celltypist", FakeCellTypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", FakeCellTypist.models)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)
    adata = ad.AnnData(
        X=np.log1p(np.array([[1, 2], [3, 4]], dtype=float)),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    model = tmp_path / "mock.pkl"
    model.write_text("mock", encoding="utf-8")
    tables = tmp_path / "tables"
    tables.mkdir()

    row, artifacts = scrna_smoke._run_celltypist_backend(
        adata=adata,
        tables=tables,
        analysis_level="validated_backend",
        model_path=model,
    )

    assert row["status"] == "failed"
    assert "no_matching_model_genes" in row["reason"]
    assert "ready" not in Path(artifacts["celltypist_annotation"]).read_text(encoding="utf-8")


def test_celltypist_backend_reports_missing_model_without_implicit_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ad = pytest.importorskip("anndata")

    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "celltypist":
            return object()
        return original_find_spec(name)

    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)
    adata = ad.AnnData(
        X=np.log1p(np.array([[1, 2], [3, 4]], dtype=float)),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    tables = tmp_path / "tables"
    tables.mkdir()

    row, artifacts = scrna_smoke._run_celltypist_backend(
        adata=adata,
        tables=tables,
        analysis_level="validated_backend",
        model_path=tmp_path / "missing_model.pkl",
    )

    assert row["status"] == "skipped"
    assert row["reason"] == "model_not_configured:celltypist_model_required_to_avoid_implicit_download"
    warning = pd.read_csv(artifacts["annotation_warning"], sep="\t")
    assert "implicit_download" in warning.loc[0, "warning"]


def test_celltypist_backend_reports_damaged_model(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ad = pytest.importorskip("anndata")

    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "celltypist":
            return object()
        return original_find_spec(name)

    class FakeCellTypist:
        class models:
            class Model:
                @staticmethod
                def load(path: str):
                    raise ValueError("not a valid CellTypist model")

    monkeypatch.setitem(sys.modules, "celltypist", FakeCellTypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", FakeCellTypist.models)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)
    model = tmp_path / "damaged.pkl"
    model.write_text("not a model", encoding="utf-8")
    tables = tmp_path / "tables"
    tables.mkdir()
    adata = ad.AnnData(
        X=np.log1p(np.array([[1, 2], [3, 4]], dtype=float)),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    row, artifacts = scrna_smoke._run_celltypist_backend(
        adata=adata,
        tables=tables,
        analysis_level="validated_backend",
        model_path=model,
    )

    assert row["status"] == "failed"
    assert "not a valid CellTypist model" in row["reason"]
    warning = pd.read_csv(artifacts["annotation_warning"], sep="\t")
    assert "not a valid CellTypist model" in warning.loc[0, "warning"]


def test_celltypist_backend_rejects_scaled_negative_expression(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ad = pytest.importorskip("anndata")

    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "celltypist":
            return object()
        return original_find_spec(name)

    class FakeCellTypist:
        class models:
            class Model:
                @staticmethod
                def load(path: str):
                    return types.SimpleNamespace(features=["G1", "G2"])

        @staticmethod
        def annotate(adata, *, model, majority_voting: bool, over_clustering=None):
            raise AssertionError("annotate should not be called for negative scaled input")

    monkeypatch.setitem(sys.modules, "celltypist", FakeCellTypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", FakeCellTypist.models)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)
    model = tmp_path / "mock.pkl"
    model.write_text("mock", encoding="utf-8")
    tables = tmp_path / "tables"
    tables.mkdir()
    adata = ad.AnnData(
        X=np.array([[1.0, -0.5], [0.2, 0.1]], dtype=float),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    row, _ = scrna_smoke._run_celltypist_backend(
        adata=adata,
        tables=tables,
        analysis_level="validated_backend",
        model_path=model,
    )

    assert row["status"] == "failed"
    assert "contains_negative_values" in row["reason"]


def test_annotate_scrna_celltypist_cli_runs_on_existing_h5ad(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ad = pytest.importorskip("anndata")

    original_find_spec = scrna_smoke.importlib.util.find_spec

    def fake_find_spec(name: str):
        if name == "celltypist":
            return object()
        return original_find_spec(name)

    class FakeCellTypist:
        __version__ = "fake-1.0"

        class models:
            class Model:
                @staticmethod
                def load(path: str):
                    return types.SimpleNamespace(features=["G1", "G2", "G3", "G4"])

        @staticmethod
        def annotate(adata, *, model, majority_voting: bool, over_clustering=None):
            labels = pd.DataFrame(
                {
                    "predicted_labels": ["T cells", "B cells", "Monocytes"],
                    "majority_voting": ["T cells", "B cells", "Monocytes"],
                },
                index=adata.obs_names,
            )
            probabilities = pd.DataFrame(
                {
                    "T cells": [0.90, 0.10, 0.10],
                    "B cells": [0.05, 0.92, 0.10],
                    "Monocytes": [0.05, 0.08, 0.93],
                },
                index=adata.obs_names,
            )
            decisions = pd.DataFrame(
                {
                    "T cells": [2.0, -2.0, -2.0],
                    "B cells": [-2.0, 2.1, -2.0],
                    "Monocytes": [-2.0, -2.0, 2.2],
                },
                index=adata.obs_names,
            )
            return types.SimpleNamespace(predicted_labels=labels, probability_matrix=probabilities, decision_matrix=decisions)

    monkeypatch.setitem(sys.modules, "celltypist", FakeCellTypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", FakeCellTypist.models)
    monkeypatch.setattr(scrna_smoke.importlib.util, "find_spec", fake_find_spec)
    input_h5ad = tmp_path / "full_gene.h5ad"
    adata = ad.AnnData(
        X=np.log1p(np.array([[4, 0, 0, 1], [0, 5, 1, 0], [0, 1, 5, 0]], dtype=float)),
        obs=pd.DataFrame({"cluster": ["0", "1", "2"]}, index=["cell1", "cell2", "cell3"]),
        var=pd.DataFrame(index=["G1", "G2", "G3", "G4"]),
    )
    adata.uns["log1p"] = {"base": None}
    adata.write_h5ad(input_h5ad)
    model = tmp_path / "mock.pkl"
    model.write_text("mock", encoding="utf-8")
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "annotate-scrna-celltypist",
            "--input-h5ad",
            str(input_h5ad),
            "--output-dir",
            str(tmp_path / "annotation"),
            "--model",
            str(model),
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["status"] == "ready"
    assert Path(payload["artifacts"]["celltypist_annotation"]).exists()
    assert Path(payload["artifacts"]["celltypist_annotation_manifest"]).exists()
    annotated = ad.read_h5ad(payload["artifacts"]["celltypist_annotated_h5ad"])
    assert "celltypist_raw_label" in annotated.obs
    assert annotated.uns["ultimate_celltypist"]["n_matched_model_features"] == 4


def test_validate_scrna_cli_refuses_demo_as_validated_backend(tmp_path: Path) -> None:
    _require_scrna_runtime()
    demo = create_demo_inputs(tmp_path / "demo_cli_guard", n_cells=32, n_genes=50, seed=14)
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "validate-scrna",
            "--input-path",
            str(demo["tenx_mtx"]),
            "--input-type",
            "10x_mtx",
            "--output-dir",
            str(tmp_path / "bad"),
            "--analysis-level",
            "validated_backend",
        ],
    )
    assert result.exit_code != 0
    assert "cannot be labeled as validated_backend" in result.output


def test_validate_scrna_cli_rejects_bad_input_type(tmp_path: Path) -> None:
    runner = CliRunner()
    path = tmp_path / "input.h5ad"
    path.write_text("not used", encoding="utf-8")
    result = runner.invoke(
        main,
        ["validate-scrna", "--input-path", str(path), "--input-type", "bad", "--output-dir", str(tmp_path / "out")],
    )
    assert result.exit_code != 0
    assert "Invalid value for '--input-type'" in result.output


def test_validate_scrna_cli_rejects_missing_input_path(tmp_path: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["validate-scrna", "--input-path", str(tmp_path / "missing.h5ad"), "--input-type", "h5ad", "--output-dir", str(tmp_path / "out")],
    )
    assert result.exit_code != 0
    assert "does not exist" in result.output


def test_unified_run_scrna_preset_uses_scrna_mvp_backend(tmp_path: Path) -> None:
    _require_scrna_runtime()
    demo = create_demo_inputs(tmp_path / "demo", n_cells=36, n_genes=45, seed=9)
    if not demo.get("h5ad"):
        pytest.skip("anndata is required to create h5ad demo input")
    config_path = tmp_path / "project.yaml"
    dump_yaml(
        {
            "project": {
                "name": "scrna_communication_unified",
                "organism": "human",
                "output_dir": str(tmp_path / "run"),
                "server_root": str(tmp_path),
                "run_mode": "interactive",
            },
            "samples": {"items": [{"sample_id": "S1", "condition": "control", "input_path": demo["h5ad"]}]},
            "modules": {
                "scrna": {
                    "enabled": True,
                    "preset": "communication",
                    "input_h5ad": demo["h5ad"],
                    "samplesheet": demo["samplesheet"],
                    "max_cells": 30,
                    "backends": {"communication": "liana,cellchat,nichenet"},
                }
            },
        },
        config_path,
    )

    manifest = run_pipeline_from_config(config_path)
    run_dir = Path(manifest["output_dir"])
    scrna = manifest["modules"][0]
    active_ids = {row["backend_id"] for row in scrna["backend_plan"]["active_backends"]}

    assert scrna["module"] == "scrna"
    assert scrna["backend_id"] == "scrna.mvp.validate_scrna"
    assert "scrna.communication.cellchat_optional" in active_ids
    assert "scrna.communication.liana" in active_ids
    assert "scrna.communication.nichenet_optional" in active_ids
    assert any(row["backend_id"] == "scrna.communication.cellchat_optional" for row in scrna["backend_execution_rows"])
    assert any(row["backend_id"] == "scrna.communication.nichenet_optional" for row in scrna["backend_execution_rows"])
    assert (run_dir / "objects" / "scrna_mvp.h5ad").exists()
    assert (run_dir / "results" / "tables" / "backend_execution_manifest.json").exists()
    assert (run_dir / "results" / "tables" / "advanced_backend_execution_manifest.json").exists()
    advanced = json.loads((run_dir / "results" / "tables" / "advanced_backend_execution_manifest.json").read_text(encoding="utf-8"))
    cellchat_rows = [row for row in advanced["rows"] if row["backend_id"] == "scrna.communication.cellchat_optional"]
    assert cellchat_rows and cellchat_rows[0]["execution_status"] != "registered_active"
    nichenet_rows = [row for row in advanced["rows"] if row["backend_id"] == "scrna.communication.nichenet_optional"]
    assert nichenet_rows and nichenet_rows[0]["execution_status"] != "registered_active"
    methods = (run_dir / "reports" / "methods.md").read_text(encoding="utf-8")
    assert "Advanced Backend Execution" in methods
    assert "scrna.communication.cellchat_optional" in methods
    assert "scrna.communication.nichenet_optional" in methods


def test_scrna_nichenet_backend_runs_with_reviewed_labels_and_resource(tmp_path: Path) -> None:
    _require_scrna_runtime()
    demo = create_demo_inputs(tmp_path / "demo_nichenet", n_cells=60, n_genes=70, seed=22)
    samples = pd.read_csv(demo["samplesheet"], sep="\t")
    labels = ["reviewed_T_like", "reviewed_B_like", "reviewed_myeloid_like"]
    samples["reviewed_cell_type"] = [labels[idx % len(labels)] for idx in range(len(samples))]
    samplesheet = tmp_path / "reviewed_labels.tsv"
    samples.to_csv(samplesheet, sep="\t", index=False)
    resource = tmp_path / "nichenet_resource.tsv"
    pd.DataFrame(
        [
            {"ligand": "IL6", "target": "MKI67", "weight": 0.9},
            {"ligand": "IL6", "target": "TOP2A", "weight": 0.8},
            {"ligand": "VEGFA", "target": "LDHA", "weight": 0.7},
            {"ligand": "VEGFA", "target": "VIM", "weight": 0.6},
        ]
    ).to_csv(resource, sep="\t", index=False)

    manifest = run_scrna_validation(
        input_path=Path(demo["tenx_mtx"]),
        input_type="10x_mtx",
        output_dir=tmp_path / "run_nichenet",
        samplesheet=samplesheet,
        max_cells=60,
        nichenet_resource=resource,
    )

    row = next(row for row in manifest["backend_status"] if row["backend_id"] == "scrna.communication.nichenet_optional")
    assert row["status"] == "ready"
    activity = pd.read_csv(tmp_path / "run_nichenet" / "results" / "tables" / "nichenet_ligand_activity.tsv", sep="\t")
    assert {"ligand", "activity_score", "warning"}.issubset(activity.columns)
    assert not activity.empty
    backend_manifest = json.loads((tmp_path / "run_nichenet" / "results" / "tables" / "nichenet_backend_manifest.json").read_text(encoding="utf-8"))
    assert backend_manifest["delivery_allowed"] is False


def test_cellchat_r_script_uses_raw_expression_for_small_validation_sets() -> None:
    source = (Path(__file__).parents[1] / "src" / "ultimate" / "scrna_smoke.py").read_text(encoding="utf-8")

    assert "computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)" in source


def test_scrna_communication_input_does_not_auto_select_velocity() -> None:
    config = {
        "modules": {
            "scrna": {
                "enabled": True,
                "preset": "communication",
                "input_path": "/shared/example/pbmc3k",
                "input_type": "10x_mtx",
                "backends": {"communication": "liana,cellchat"},
            }
        }
    }

    assert has_scrna_velocity_backend_config(config) is False


def test_scrna_trajectory_input_selects_velocity() -> None:
    config = {
        "modules": {
            "scrna": {
                "enabled": True,
                "preset": "trajectory",
                "input_h5ad": "/shared/example/velocity.h5ad",
            }
        }
    }

    assert has_scrna_velocity_backend_config(config) is True


def test_scrna_mvp_slurm_uses_explicit_celltypist_reference_cache() -> None:
    script = (Path(__file__).parents[1] / "slurm" / "scrna_mvp_validation.sbatch").read_text(encoding="utf-8")

    assert "CELLTYPIST_FOLDER" in script
    assert "references/celltypist" in script
    assert "Immune_All_Low.pkl" in script
    assert "--celltypist-model" in script
    assert "export CELLTYPIST_FOLDER" in script


def test_scrna_backend_rows_record_slurm_context(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _require_scrna_runtime()
    monkeypatch.setenv("SLURM_JOB_ID", "pytest-123")
    monkeypatch.setenv("SLURM_JOB_NAME", "pytest_scrna")
    demo = create_demo_inputs(tmp_path / "demo_slurm_context", n_cells=32, n_genes=50, seed=15)

    manifest = run_scrna_validation(
        input_path=Path(demo["h5ad"]),
        input_type="h5ad",
        output_dir=tmp_path / "run_slurm_context",
        samplesheet=Path(demo["samplesheet"]),
        max_cells=32,
    )

    for row in manifest["backend_status"]:
        assert row["backend_slurm_job_id"] == "pytest-123"
        assert row["backend_slurm_job_name"] == "pytest_scrna"


def test_pseudobulk_de_backend_executes_when_r_backend_available(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    tables = tmp_path / "results" / "tables"
    tables.mkdir(parents=True)
    counts = pd.DataFrame(
        {
            "feature_id": ["GeneA", "GeneB"],
            "S1__control__cluster0": [10, 20],
            "S2__control__cluster0": [12, 18],
            "S3__case__cluster0": [30, 5],
            "S4__case__cluster0": [28, 6],
        }
    )
    design = pd.DataFrame(
        {
            "pseudobulk_id": ["S1__control__cluster0", "S2__control__cluster0", "S3__case__cluster0", "S4__case__cluster0"],
            "sample_id": ["S1", "S2", "S3", "S4"],
            "condition": ["control", "control", "case", "case"],
            "cluster": ["0", "0", "0", "0"],
        }
    )
    counts.to_csv(tables / "pseudobulk_counts.tsv", sep="\t", index=False)
    design.to_csv(tables / "pseudobulk_design.tsv", sep="\t", index=False)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_rscript = fake_bin / "Rscript"
    fake_rscript.write_text(
        """#!/usr/bin/env python3
import json
import sys
from pathlib import Path
if "-e" in sys.argv:
    print("DESeq2\\nedgeR\\njsonlite")
    raise SystemExit(0)
args = sys.argv[1:]
tables = Path(args[args.index("--tables-dir") + 1])
figures = tables.parent / "figures"
objects = tables.parent.parent / "objects"
figures.mkdir(parents=True, exist_ok=True)
objects.mkdir(parents=True, exist_ok=True)
(tables / "pseudobulk_de_results.tsv").write_text("feature_id\\tcluster\\tlog2FoldChange\\tpvalue\\tpadj\\tbackend_id\\tbackend_method\\tanalysis_level\\nGeneA\\t0\\t1.5\\t0.01\\t0.05\\tscrna.pseudobulk.deseq2_edger\\tDESeq2\\tvalidated_backend\\n")
(tables / "pseudobulk_de_backend_status.tsv").write_text("backend_id\\tcluster\\tstatus\\tanalysis_level\\tbackend_method\\nscrna.pseudobulk.deseq2_edger\\t0\\tready\\tvalidated_backend\\tDESeq2\\n")
(tables / "pseudobulk_de_backend_versions.tsv").write_text("package\\tversion\\nDESeq2\\t1.0\\nedgeR\\t1.0\\njsonlite\\t1.0\\n")
(tables / "pseudobulk_de_backend_manifest.json").write_text(json.dumps({"backend_id":"scrna.pseudobulk.deseq2_edger","status":"ready","analysis_level":"validated_backend"}))
(figures / "pseudobulk_de_volcano.png").write_text("png")
(objects / "scrna_pseudobulk_de_backend.rds").write_text("rds")
""",
        encoding="utf-8",
    )
    fake_rscript.chmod(0o755)
    import os

    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}")

    row, artifacts = scrna_smoke._run_pseudobulk_de_backend(tables=tables, analysis_level="validated_backend")

    assert row["status"] == "ready"
    assert Path(artifacts["pseudobulk_de_backend_manifest"]).exists()
    assert Path(artifacts["pseudobulk_de_volcano"]).exists()
    assert Path(artifacts["pseudobulk_de_rds"]).exists()


def test_validate_scrna_pseudobulk_de_cli_invokes_paired_backend(tmp_path: Path) -> None:
    fake_rscript = tmp_path / "Rscript"
    fake_rscript.write_text(
        """#!/usr/bin/env python3
import json
import sys
from pathlib import Path
args = sys.argv[1:]
tables = Path(args[args.index("--tables-dir") + 1])
figures = tables.parent / "figures"
objects = tables.parent.parent / "objects"
tables.mkdir(parents=True, exist_ok=True)
figures.mkdir(parents=True, exist_ok=True)
objects.mkdir(parents=True, exist_ok=True)
(tables / "pseudobulk_de_results.tsv").write_text("feature_id\\tlog2FoldChange\\tpvalue\\tpadj\\tcell_type\\tbackend_id\\tbackend_method\\nGeneA\\t1.0\\t0.01\\t0.05\\tCD4 T cells\\tscrna.pseudobulk.edger_paired\\tedgeR\\n")
(tables / "pseudobulk_de_backend_status.tsv").write_text("backend_id\\tcell_type\\tstatus\\tbackend_method\\n")
(tables / "pseudobulk_de_backend_versions.tsv").write_text("package\\tversion\\nedgeR\\t4.4.0\\n")
(figures / "pseudobulk_de_volcano.png").write_text("png")
(objects / "scrna_pseudobulk_de_backend.rds").write_text("rds")
(tables / "pseudobulk_de_backend_manifest.json").write_text(json.dumps({
    "backend_id": "scrna.pseudobulk.edger_paired",
    "status": "ready",
    "backend_method": "edgeR",
    "target_cell_type": "CD4 T cells",
}))
""",
        encoding="utf-8",
    )
    fake_rscript.chmod(0o755)
    sce = tmp_path / "dummy.Rda"
    sce.write_text("not used by fake backend", encoding="utf-8")

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "validate-scrna-pseudobulk-de",
            "--sce-rda",
            str(sce),
            "--output-dir",
            str(tmp_path / "run"),
            "--rscript",
            str(fake_rscript),
            "--cell-type",
            "CD4 T cells",
            "--control",
            "ctrl",
            "--case",
            "stim",
            "--multiplet-column",
            "multiplets",
            "--multiplet-keep",
            "singlet",
            "--min-cells",
            "50",
            "--min-pairs",
            "4",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["backend_id"] == "scrna.pseudobulk.edger_paired"
    assert payload["status"] == "ready"
    assert payload["target_cell_type"] == "CD4 T cells"
    assert Path(payload["ultimate_cli_log"]).exists()


def test_scrna_paired_pseudobulk_edger_r_regression() -> None:
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript is required for paired edgeR pseudobulk regression")
    package_check = subprocess.run(
        [
            rscript,
            "-e",
            "pkgs <- c('edgeR','limma','SingleCellExperiment','Matrix','jsonlite'); ok <- vapply(pkgs, requireNamespace, logical(1), quietly=TRUE); quit(status=ifelse(all(ok), 0, 1))",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    if package_check.returncode != 0:
        pytest.skip("edgeR/SingleCellExperiment R packages are required for paired pseudobulk regression")
    script = Path(__file__).resolve().parent / "R" / "test_scrna_paired_pseudobulk_edger.R"
    completed = subprocess.run([rscript, str(script)], text=True, capture_output=True, check=False)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "scrna_paired_pseudobulk_edger_regression_pass" in completed.stdout


def _require_scrna_runtime() -> None:
    pytest.importorskip("scanpy")
    pytest.importorskip("anndata")
    pytest.importorskip("h5py")
    pytest.importorskip("scipy")


def _assert_mvp_outputs(manifest: dict) -> None:
    run_dir = Path(manifest["output_dir"])
    assert manifest["status"] == "ready"
    assert manifest["objects"]["h5ad"].endswith("objects/scrna_mvp.h5ad")
    assert Path(manifest["objects"]["h5ad"]).exists()
    assert not (run_dir / "objects" / "scrna_smoke_validated.h5ad").exists()
    assert manifest["cell_type_annotation_status"] == "placeholder_not_cell_type"
    assert manifest["pseudobulk_de_status"] == "design_ready_matrix_only"
    ad = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    adata = ad.read_h5ad(run_dir / "objects" / "scrna_mvp.h5ad")
    assert "counts" in adata.layers
    assert "full_gene_h5ad" in manifest["objects"]
    assert Path(manifest["objects"]["full_gene_h5ad"]).exists()
    assert Path(manifest["gene_retention_summary"]).exists()
    assert Path(manifest["matrix_semantics"]).exists()
    assert Path(manifest["hvg_genes"]).exists()
    counts = adata.layers["counts"]
    count_values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    assert count_values.size and np.nanmin(count_values) >= 0
    assert "backend_execution_manifest" in manifest
    assert "backend_execution_status" in manifest
    assert Path(manifest["backend_execution_manifest"]).exists()
    handoff_text = (run_dir / "results/tables/pseudobulk_deseq2_edgeR_handoff.R").read_text(encoding="utf-8")
    assert "Guarded backend" in handoff_text
    assert "not a formal DE result" in handoff_text
    assert "TODO" not in handoff_text
    backend_ids = {row["backend_id"] for row in manifest["backend_status"]}
    assert {
        "scrna.qc.scrublet",
        "scrna.annotation.celltypist",
        "scrna.functional.decoupler_gseapy",
        "scrna.communication.liana",
        "scrna.communication.cellchat_optional",
        "scrna.communication.nichenet_optional",
        "scrna.pseudobulk.deseq2_edger",
    }.issubset(backend_ids)
    assert Path(manifest["raw_qc_manifest"]).exists()
    for relative in [
        "results/tables/qc_metrics.tsv",
        "results/tables/marker_genes.tsv",
        "results/tables/de_condition.tsv",
        "results/tables/cell_type_composition.tsv",
        "results/tables/basic_enrichment.tsv",
        "results/tables/cell_type_annotation_placeholder.tsv",
        "results/tables/pseudobulk_counts.tsv",
        "results/tables/pseudobulk_design.tsv",
        "results/tables/pseudobulk_feature_metadata.tsv",
        "results/tables/backend_execution.tsv",
        "results/tables/backend_execution_manifest.json",
        "results/tables/gene_retention_summary.tsv",
        "results/tables/hvg_genes.tsv",
        "results/tables/matrix_semantics.tsv",
        "results/tables/doublet_scores.tsv",
        "results/tables/doublet_summary.tsv",
        "results/tables/celltypist_annotation.tsv",
        "results/tables/annotation_confidence.tsv",
        "results/tables/annotation_warning.tsv",
        "results/tables/signature_scores.tsv",
        "results/tables/pathway_activity.tsv",
        "results/tables/tf_activity.tsv",
        "results/tables/functional_backend_status.tsv",
        "results/tables/liana_interactions.tsv",
        "results/tables/communication_network.tsv",
        "results/tables/communication_backend_status.tsv",
        "results/tables/cellchat_interactions.tsv",
        "results/tables/cellchat_pathway_summary.tsv",
        "results/tables/cellchat_backend_status.tsv",
        "results/tables/cellchat_backend_manifest.json",
        "results/tables/nichenet_ligand_activity.tsv",
        "results/tables/nichenet_ligand_target_links.tsv",
        "results/tables/nichenet_backend_status.tsv",
        "results/tables/nichenet_backend_manifest.json",
        "results/tables/pseudobulk_de_backend_status.tsv",
        "results/tables/pseudobulk_de_results.tsv",
        "results/tables/pseudobulk_deseq2_edgeR_handoff.R",
        "results/figures/qc_violin.png",
        "results/figures/pca_condition.png",
        "results/figures/umap_cluster_condition.png",
        "results/figures/volcano_condition.png",
        "results/figures/communication_dotplot.png",
        "results/figures/cellchat_network.png",
        "results/figures/nichenet_activity.png",
        "reports/report.md",
        "reports/report.html",
        "run_manifest.json",
    ]:
        path = run_dir / relative
        assert path.exists(), relative
        assert path.stat().st_size > 0, relative
    report = (run_dir / "reports" / "report.md").read_text(encoding="utf-8")
    assert "analysis_level" in report
    assert "delivery_allowed" in report
    assert "backend 执行摘要" in report
    assert "cluster placeholder" in report
    assert "LIANA" in report
    manifest_from_disk = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest_from_disk["analysis_level"] == manifest["analysis_level"]
    annotation = Path(run_dir / "results/tables/cell_type_annotation_placeholder.tsv").read_text(encoding="utf-8")
    assert "placeholder_not_cell_type" in annotation
    design = Path(run_dir / "results/tables/pseudobulk_design.tsv").read_text(encoding="utf-8")
    assert "analysis_level" in design
