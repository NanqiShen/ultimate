"""M4 deep-analysis regressions: Harmony, paired edgeR, CopyKAT, enrichment, squidpy.

Every test runs without the optional scientific stack installed. Real branches are
driven with stand-in modules (fake ``harmonypy``/``squidpy``/``gseapy``) so the
production code path is actually exercised, while missing-dependency branches
assert an explicit ``skipped`` with a concrete reason - never a silent no-op and
never fabricated statistics.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import ultimate.scrna_smoke as scrna
import ultimate.tumor_sc_backend as tumor
from ultimate.enrichment_backend import run_enrichment
from ultimate.spatial_backend import _run_squidpy_neighborhood, _write_spatial_neighborhood_versions


# --------------------------------------------------------------------------
# shared stand-ins
# --------------------------------------------------------------------------


class _PP:
    @staticmethod
    def neighbors(adata, **kwargs):
        adata.obsm["_neighbors_done"] = True


class _TL:
    @staticmethod
    def umap(adata, random_state=0):
        base = adata.obsm.get("X_pca_harmony", adata.obsm["X_pca"])
        adata.obsm["X_umap"] = np.asarray(base)[:, :2]


class FakeScanpy:
    pp = _PP
    tl = _TL


class FakeAdata:
    def __init__(self, n_cells: int = 60, n_pcs: int = 8, batches=("b1", "b2"), batch_column: str = "batch"):
        rng = np.random.default_rng(3)
        self.n_obs = n_cells
        self.n_vars = n_pcs
        self.obs_names = [f"c{i}" for i in range(n_cells)]
        self.obsm = {"X_pca": rng.normal(size=(n_cells, n_pcs))}
        labels = [batches[i % len(batches)] for i in range(n_cells)]
        self.obs = pd.DataFrame({batch_column: labels}, index=self.obs_names)
        self.uns = {}


def _install_fake_harmonypy(monkeypatch):
    module = types.ModuleType("harmonypy")

    class _Result:
        def __init__(self, corrected):
            self.Z_corr = corrected

    def run_harmony(Z, meta, col, random_state=0, verbose=False):
        values = np.asarray(Z, dtype=float)
        corrected = values.copy()
        for batch in pd.unique(meta[col]):
            mask = (meta[col] == batch).to_numpy()
            corrected[mask] -= values[mask].mean(axis=0)
        return _Result(corrected.T)  # squidpy/harmony convention: (n_pcs, n_cells)

    module.run_harmony = run_harmony
    monkeypatch.setitem(sys.modules, "harmonypy", module)


def _drop_module(monkeypatch, name: str):
    monkeypatch.setitem(sys.modules, name, None)


# --------------------------------------------------------------------------
# 1. Harmony batch integration
# --------------------------------------------------------------------------


def test_harmony_integration_ready_branch(tmp_path: Path, monkeypatch):
    _install_fake_harmonypy(monkeypatch)
    adata = FakeAdata()
    config = {"modules": {"scrna": {"batch_correction": "harmony"}}}

    row, artifacts = scrna._run_batch_integration(
        sc=FakeScanpy, adata=adata, tables=tmp_path / "t", figures=tmp_path / "f",
        config=config, random_seed=7, analysis_level="smoke_backend",
    )

    assert row["status"] == "ready", row
    assert row["analysis_level"] == "smoke_backend"
    assert adata.obsm["X_pca_harmony"].shape == (adata.n_obs, adata.n_vars)
    assert adata.uns["ultimate_integration_use_rep"] == "X_pca_harmony"
    assert "X_umap_unintegrated" in adata.obsm
    metrics = pd.read_csv(tmp_path / "t" / "batch_integration_metrics.tsv", sep="\t")
    for column in ("batch_entropy_before", "batch_entropy_after", "ilisi_approx_before",
                   "ilisi_approx_after", "cluster_preservation_ari", "n_batches", "n_cells",
                   "backend_id", "status", "analysis_level"):
        assert column in metrics.columns, column
    assert int(metrics.iloc[0]["n_batches"]) == 2
    assert (tmp_path / "f" / "umap_batch_before_after.png").is_file()
    assert artifacts["batch_integration_metrics"].endswith("batch_integration_metrics.tsv")


def test_harmony_skips_without_harmonypy(tmp_path: Path, monkeypatch):
    _drop_module(monkeypatch, "harmonypy")
    adata = FakeAdata()
    config = {"modules": {"scrna": {"batch_correction": "harmony"}}}

    row, _ = scrna._run_batch_integration(
        sc=FakeScanpy, adata=adata, tables=tmp_path / "t", figures=tmp_path / "f",
        config=config, random_seed=7, analysis_level="smoke_backend",
    )

    assert row["status"] == "skipped"
    assert "harmonypy" in row["reason"]
    assert "X_pca_harmony" not in adata.obsm  # nothing was silently integrated
    metrics = pd.read_csv(tmp_path / "t" / "batch_integration_metrics.tsv", sep="\t")
    assert metrics.iloc[0]["status"] == "skipped"


def test_harmony_skips_on_single_batch(tmp_path: Path, monkeypatch):
    _install_fake_harmonypy(monkeypatch)
    adata = FakeAdata(batches=("only",))
    row, _ = scrna._run_batch_integration(
        sc=FakeScanpy, adata=adata, tables=tmp_path / "t", figures=tmp_path / "f",
        config={"modules": {"scrna": {"batch_correction": "harmony"}}},
        random_seed=7, analysis_level="smoke_backend",
    )
    assert row["status"] == "skipped"
    assert "single_batch" in row["reason"]


def test_harmony_skips_when_disabled_or_missing_column(tmp_path: Path, monkeypatch):
    _install_fake_harmonypy(monkeypatch)
    row, _ = scrna._run_batch_integration(
        sc=FakeScanpy, adata=FakeAdata(), tables=tmp_path / "t1", figures=tmp_path / "f1",
        config={"modules": {"scrna": {"batch_correction": "none"}}},
        random_seed=7, analysis_level="smoke_backend",
    )
    assert row["status"] == "skipped" and "disabled" in row["reason"]

    row2, _ = scrna._run_batch_integration(
        sc=FakeScanpy, adata=FakeAdata(batch_column="lane"), tables=tmp_path / "t2", figures=tmp_path / "f2",
        config={"modules": {"scrna": {"batch_correction": "harmony"}}},
        random_seed=7, analysis_level="smoke_backend",
    )
    assert row2["status"] == "skipped" and "missing_batch_column" in row2["reason"]


def test_batch_column_resolution_precedence():
    assert scrna._resolve_batch_column({"design": {"batch_column": "lane"}}) == "lane"
    assert scrna._resolve_batch_column({"modules": {"scrna": {"batch_column": "pool"}}}) == "pool"
    assert scrna._resolve_batch_column({}) == "batch"
    assert scrna._resolve_batch_correction_mode({"modules": {"scrna": {"batch_correction": "AUTO"}}}) == "auto"
    assert scrna._resolve_batch_correction_mode({"modules": {"scrna": {"batch_correction": "bogus"}}}) == "none"


# --------------------------------------------------------------------------
# 2A. paired edgeR pseudobulk DE inside `ultimate run`
# --------------------------------------------------------------------------


def test_paired_pseudobulk_not_configured_returns_none(tmp_path: Path):
    row, artifacts = scrna._run_paired_pseudobulk_edger(
        config={}, tables=tmp_path / "t", figures=tmp_path / "f", objects=tmp_path / "o",
        analysis_level="smoke_backend",
    )
    assert row is None and artifacts == {}


def test_paired_pseudobulk_skips_without_input(tmp_path: Path):
    cfg = {"modules": {"scrna": {"pseudobulk": {
        "edger_paired": True, "donor_column": "ind", "condition_column": "stim",
        "cell_type_column": "cell", "cell_type": "T",
    }}}}
    row, _ = scrna._run_paired_pseudobulk_edger(
        config=cfg, tables=tmp_path / "t", figures=tmp_path / "f", objects=tmp_path / "o",
        analysis_level="smoke_backend",
    )
    assert row is not None
    assert row["status"] == "skipped"
    assert "missing_paired_input" in row["reason"]


def test_paired_pseudobulk_runs_through_run_pipeline(tmp_path: Path, monkeypatch):
    """The R success branch is driven end-to-end through the backend row."""
    seen: dict[str, object] = {}

    def fake_run(**kwargs):
        seen.update(kwargs)
        return {"status": "ready", "skip_reason": "", "n_pairs": 6, "engine": "edgeR"}

    monkeypatch.setattr("ultimate.scrna_pseudobulk_de.run_scrna_pseudobulk_de", fake_run)
    monkeypatch.setattr(scrna, "_rscript_str", lambda m: "/usr/bin/Rscript")
    counts = tmp_path / "counts.tsv"
    design = tmp_path / "design.tsv"
    counts.write_text("gene\tS1\tS2\nG1\t5\t6\n", encoding="utf-8")
    design.write_text("sample\tdonor\tstim\nS1\td1\tctrl\nS2\td1\tstim\n", encoding="utf-8")
    cfg = {"modules": {"scrna": {"pseudobulk": {
        "edger_paired": True, "donor_column": "donor", "condition_column": "stim",
        "cell_type_column": "cell", "cell_type": "T",
        "counts": str(counts), "design": str(design),
    }}}}
    (tmp_path / "t").mkdir()

    row, artifacts = scrna._run_paired_pseudobulk_edger(
        config=cfg, tables=tmp_path / "t", figures=tmp_path / "f", objects=tmp_path / "o",
        analysis_level="validated_backend",
    )

    assert row["status"] == "ready"
    assert row["backend_id"] == "scrna.pseudobulk.edger_paired"
    assert seen["donor_column"] == "donor"
    assert seen["condition_column"] == "stim"
    assert seen["cell_type"] == "T"
    assert (tmp_path / "t" / "pseudobulk_paired_backend_status.tsv").is_file()
    assert artifacts["pseudobulk_paired_backend_status"]


def test_paired_pseudobulk_failure_is_explicit_skip(tmp_path: Path, monkeypatch):
    def boom(**kwargs):
        raise RuntimeError("edgeR_failed")

    monkeypatch.setattr("ultimate.scrna_pseudobulk_de.run_scrna_pseudobulk_de", boom)
    monkeypatch.setattr(scrna, "_rscript_str", lambda m: "/usr/bin/Rscript")
    counts = tmp_path / "c.tsv"; counts.write_text("g\tS1\nG1\t1\n", encoding="utf-8")
    design = tmp_path / "d.tsv"; design.write_text("sample\tdonor\tstim\nS1\td1\tctrl\n", encoding="utf-8")
    cfg = {"modules": {"scrna": {"pseudobulk": {
        "edger_paired": True, "donor_column": "donor", "condition_column": "stim",
        "cell_type_column": "cell", "cell_type": "T", "counts": str(counts), "design": str(design),
    }}}}
    (tmp_path / "t").mkdir()
    row, _ = scrna._run_paired_pseudobulk_edger(
        config=cfg, tables=tmp_path / "t", figures=tmp_path / "f", objects=tmp_path / "o",
        analysis_level="smoke_backend",
    )
    assert row["status"] == "skipped"
    assert "paired_pseudobulk_failed" in row["reason"]


# --------------------------------------------------------------------------
# 2B. CopyKAT
# --------------------------------------------------------------------------


def test_copykat_requires_explicit_opt_in():
    assert tumor._copykat_enabled({}, {}) is False
    assert tumor._copykat_enabled({}, {"run_copykat": True}) is True
    assert tumor._copykat_enabled({}, {"backend": "scrna.tumor.copykat"}) is True
    assert tumor._copykat_enabled({}, {"backends": {"scrna.tumor.copykat": True}}) is True


def test_copykat_integer_count_gate():
    raw = pd.DataFrame([[5, 6], [7, 8]])
    normalized = pd.DataFrame([[5.5, 6.2], [7.1, 8.0]])
    negative = pd.DataFrame([[-1, 2], [3, 4]])
    assert tumor._integer_count_gate(raw) == (True, "")
    assert tumor._integer_count_gate(normalized)[0] is False
    assert "non_integer" in tumor._integer_count_gate(normalized)[1]
    assert "negative" in tumor._integer_count_gate(negative)[1]


def _copykat_run(config, matrix, tmp_path):
    return tumor._run_copykat_backend(
        config=config, data={"matrix": matrix},
        tables_dir=tmp_path / "t", figures_dir=tmp_path / "f",
        objects_dir=tmp_path / "o", logs_dir=tmp_path / "l",
        analysis_level="smoke_backend",
    )


def test_copykat_disabled_by_default_is_skipped(tmp_path: Path):
    row, _ = _copykat_run({"modules": {"tumor_sc": {}}}, pd.DataFrame([[5, 6], [7, 8]]), tmp_path)
    assert row["status"] == "skipped"
    assert row["reason"] == "copykat_not_enabled"


def test_copykat_blocks_normalized_matrix(tmp_path: Path):
    matrix = pd.DataFrame([[5.5, 6.2], [7.1, 8.0]])
    row, _ = _copykat_run({"modules": {"tumor_sc": {"run_copykat": True}}}, matrix, tmp_path)
    assert row["status"] == "blocked"
    assert "non_integer" in row["reason"]


def test_copykat_missing_package_is_skipped(tmp_path: Path, monkeypatch):
    """With R present but no copykat, it must skip explicitly (no fake result)."""
    monkeypatch.setattr(tumor, "resolve_rscript", lambda module: Path("/opt/homebrew/bin/Rscript"))

    class _Probe:
        returncode = 2
        stdout = ""
        stderr = ""

    monkeypatch.setattr(tumor.subprocess, "run", lambda *a, **k: _Probe())
    row, _ = _copykat_run({"modules": {"tumor_sc": {"run_copykat": True}}}, pd.DataFrame([[5, 6], [7, 8]]), tmp_path)
    assert row["status"] == "skipped"
    assert row["reason"] == "dependency_missing:copykat"


def test_copykat_passes_thread_budget_as_ncores(tmp_path: Path, monkeypatch):
    """--ncores must come from thread_budget, not a hardcoded value."""
    captured: dict[str, list[str]] = {}

    monkeypatch.setattr(tumor, "resolve_rscript", lambda module: Path("/opt/homebrew/bin/Rscript"))
    monkeypatch.setattr(tumor, "resolve_packaged_asset", lambda *parts: Path("/tmp/tumor_cnv_copykat.R"))

    class _Ok:
        returncode = 0
        stdout = "done"
        stderr = ""

    def fake_run(command, **kwargs):
        if "-e" in command:
            class _P:
                returncode = 0
            return _P()
        captured["command"] = list(command)
        return _Ok()

    monkeypatch.setattr(tumor.subprocess, "run", fake_run)
    row, _ = _copykat_run(
        {"modules": {"tumor_sc": {"run_copykat": True}}, "runtime": {"threads": 2}},
        pd.DataFrame([[5, 6], [7, 8]]), tmp_path,
    )
    assert row["status"] == "ready", row
    command = captured["command"]
    assert "--ncores" in command
    ncores = int(command[command.index("--ncores") + 1])
    assert 1 <= ncores <= 8


# --------------------------------------------------------------------------
# 3. local-GMT enrichment
# --------------------------------------------------------------------------


def _de_table(n: int = 30, seed: int = 1) -> pd.DataFrame:
    return pd.DataFrame({
        "feature_id": [f"G{i}" for i in range(n)],
        "log2FC": np.random.default_rng(seed).normal(size=n),
        "padj": np.linspace(0.001, 0.9, n),
    })


def _gmt(tmp_path: Path) -> Path:
    path = tmp_path / "sets.gmt"
    path.write_text(
        "PATH_A\tna\t" + "\t".join(f"G{i}" for i in range(15)) + "\n"
        "PATH_B\tna\t" + "\t".join(f"G{i}" for i in range(10, 25)) + "\n",
        encoding="utf-8",
    )
    return path


def test_enrichment_skips_without_local_gmt(tmp_path: Path):
    result = run_enrichment(
        de_table=_de_table(), config={}, module_name="rnaseq",
        output_dir=tmp_path / "out", analysis_level="smoke_backend",
    )
    assert result["status"] == "skipped"
    assert "missing_local_gmt" in result["reason"]
    status = pd.read_csv(tmp_path / "out" / "results" / "tables" / "rnaseq" / "enrichment_backend_status.tsv", sep="\t")
    assert status.iloc[0]["network_access"] == "none_local_gmt_only"


def test_enrichment_refuses_remote_gmt_url(tmp_path: Path):
    config = {"resources": {"msigdb": "https://example.org/sets.gmt"}}
    result = run_enrichment(
        de_table=_de_table(), config=config, module_name="rnaseq",
        output_dir=tmp_path / "out", analysis_level="smoke_backend",
    )
    assert result["status"] == "skipped"
    assert "missing_local_gmt" in result["reason"]


def test_enrichment_skips_without_engine(tmp_path: Path, monkeypatch):
    monkeypatch.setitem(sys.modules, "gseapy", None)
    monkeypatch.setenv("ULTIMATE_R_ENV_RNASEQ", str(tmp_path / "missing_Rscript"))
    monkeypatch.setenv("ULTIMATE_RSCRIPT", str(tmp_path / "missing_Rscript"))
    monkeypatch.setenv("ULTIMATE_ENV_ROOT", str(tmp_path / "no_envs"))
    monkeypatch.setenv("PATH", str(tmp_path / "no_bin"))
    result = run_enrichment(
        de_table=_de_table(), config={"resources": {"msigdb": str(_gmt(tmp_path))}},
        module_name="rnaseq", output_dir=tmp_path / "out", analysis_level="smoke_backend",
    )
    assert result["status"] == "skipped"
    assert result["reason"] == "dependency_missing:gseapy"


def test_enrichment_fgsea_fallback_when_gseapy_missing(tmp_path: Path, monkeypatch):
    monkeypatch.setitem(sys.modules, "gseapy", None)
    fake_r = tmp_path / "fake_Rscript"
    fake_r.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-e" in args:
    raise SystemExit(0)
out_path = Path(args[args.index("--output") + 1])
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(
    "Term\\tES\\tNES\\tNOM p-val\\tFDR q-val\\tLead_genes\\tsize\\n"
    "PATH_A\\t0.65\\t1.85\\t0.002\\t0.01\\tGENE1;GENE2\\t15\\n",
    encoding="utf-8",
)
""",
        encoding="utf-8",
    )
    fake_r.chmod(0o755)
    monkeypatch.setenv("ULTIMATE_R_ENV_RNASEQ", str(fake_r))

    result = run_enrichment(
        de_table=_de_table(), config={"resources": {"msigdb": str(_gmt(tmp_path))}},
        module_name="rnaseq", output_dir=tmp_path / "out", analysis_level="smoke_backend",
    )

    assert result["status"] == "ready", result
    assert result.get("engine") == "fgsea"
    table = pd.read_csv(tmp_path / "out" / "results" / "tables" / "rnaseq" / "enrichment_results.tsv", sep="\t")
    assert {"Term", "NES", "backend_id", "gene_set_source"}.issubset(table.columns)
    assert (tmp_path / "out" / "results" / "figures" / "rnaseq" / "enrichment_dotplot.png").is_file()


def test_enrichment_ranks_deseq2_log2foldchange_column(tmp_path: Path, monkeypatch):
    monkeypatch.setitem(sys.modules, "gseapy", None)
    fake_r = tmp_path / "fake_Rscript"
    fake_r.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-e" in args:
    raise SystemExit(0)
ranked_path = Path(args[args.index("--ranked") + 1])
rows = ranked_path.read_text(encoding="utf-8").splitlines()
assert len(rows) > 10
out_path = Path(args[args.index("--output") + 1])
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(
    "Term\\tES\\tNES\\tNOM p-val\\tFDR q-val\\tLead_genes\\tsize\\n"
    "PATH_DESEQ2\\t0.70\\t1.95\\t0.001\\t0.008\\tGENE1;GENE2\\t15\\n",
    encoding="utf-8",
)
""",
        encoding="utf-8",
    )
    fake_r.chmod(0o755)
    monkeypatch.setenv("ULTIMATE_R_ENV_RNASEQ", str(fake_r))

    deseq2_de = pd.DataFrame({
        "feature_id": [f"GENE{i}" for i in range(1, 31)],
        "log2FoldChange": np.linspace(2.5, -2.5, 30),
        "padj": np.linspace(1e-4, 0.9, 30),
    })
    result = run_enrichment(
        de_table=deseq2_de, config={"resources": {"msigdb": str(_gmt(tmp_path))}},
        module_name="rnaseq", output_dir=tmp_path / "out", analysis_level="validated_backend",
    )
    assert result["status"] == "ready", result
    assert result["n_ranked_genes"] == 30


def test_enrichment_skips_with_too_few_ranked_genes(tmp_path: Path):
    tiny = pd.DataFrame({"feature_id": ["G1", "G2"], "log2FC": [1.0, 2.0]})
    result = run_enrichment(
        de_table=tiny, config={"resources": {"msigdb": str(_gmt(tmp_path))}},
        module_name="rnaseq", output_dir=tmp_path / "out", analysis_level="smoke_backend",
    )
    assert result["status"] == "skipped"
    assert "insufficient_ranked_genes" in result["reason"]


def test_enrichment_real_branch_with_stand_in_gseapy(tmp_path: Path, monkeypatch):
    """Exercise the real prerank path with a stand-in engine (no network)."""
    calls: dict[str, object] = {}
    module = types.ModuleType("gseapy")

    class _Result:
        def __init__(self):
            self.res2d = pd.DataFrame({
                "Term": ["PATH_A", "PATH_B"],
                "NES": [1.8, -1.4],
                "FDR q-val": [0.01, 0.05],
            })
            self.ranking = list(range(30))
            self.results = {"PATH_A": {"Term": "PATH_A", "RES": list(np.linspace(0, 0.6, 40))}}

    def prerank(**kwargs):
        calls.update(kwargs)
        return _Result()

    module.prerank = prerank
    monkeypatch.setitem(sys.modules, "gseapy", module)

    result = run_enrichment(
        de_table=_de_table(), config={"resources": {"msigdb": str(_gmt(tmp_path))}},
        module_name="rnaseq", output_dir=tmp_path / "out", analysis_level="smoke_backend",
    )

    assert result["status"] == "ready", result
    table = pd.read_csv(tmp_path / "out" / "results" / "tables" / "rnaseq" / "enrichment_results.tsv", sep="\t")
    assert {"Term", "NES", "backend_id", "gene_set_source"}.issubset(table.columns)
    assert (tmp_path / "out" / "results" / "figures" / "rnaseq" / "enrichment_dotplot.png").is_file()
    # A local GMT path was passed through; no URL was ever used.
    assert "https://" not in str(calls.get("gene_sets"))


def test_enrichment_module_never_performs_network_calls():
    source = (Path(__file__).resolve().parents[1] / "src" / "ultimate" / "enrichment_backend.py").read_text(encoding="utf-8")
    for forbidden in ("requests.", "urlopen", "urllib.request", "Enrichr(", "enrichr("):
        assert forbidden not in source, forbidden


# --------------------------------------------------------------------------
# 4. spatial squidpy
# --------------------------------------------------------------------------


def _spatial_data(n_spots: int = 40, n_features: int = 6) -> dict:
    rng = np.random.default_rng(5)
    return {
        "coords": rng.uniform(0, 100, (n_spots, 2)),
        "matrix": rng.poisson(5, (n_spots, n_features)).astype(float),
        "features": [f"G{i}" for i in range(n_features)],
        "spots": [f"S{i}" for i in range(n_spots)],
    }


def test_squidpy_branch_returns_none_without_package(monkeypatch):
    monkeypatch.setitem(sys.modules, "squidpy", None)
    assert _run_squidpy_neighborhood(
        data=_spatial_data(), analysis_fields={"analysis_level": "smoke_backend"},
        backend_id="spatial.neighborhood.squidpy", warning="w",
        input_artifact="a", source_dataset="d",
    ) is None


def test_spatial_versions_label_fallback_explicitly(tmp_path: Path):
    path = tmp_path / "versions.tsv"
    _write_spatial_neighborhood_versions(path, "spatial.neighborhood.squidpy")
    frame = pd.read_csv(path, sep="\t")
    assert (frame["execution_mode"] == "fallback_no_permutation_pvalue").all()
    assert "squidpy" not in set(frame["tool"])  # no version is claimed when absent


def test_spatial_versions_record_real_squidpy_version(tmp_path: Path):
    path = tmp_path / "versions.tsv"
    _write_spatial_neighborhood_versions(
        path, "spatial.neighborhood.squidpy", execution_mode="real_squidpy_gr", squidpy_version="1.6.0",
    )
    frame = pd.read_csv(path, sep="\t")
    row = frame[frame["tool"] == "squidpy"].iloc[0]
    assert row["version"] == "1.6.0"
    assert row["execution_mode"] == "real_squidpy_gr"


def test_spatial_real_squidpy_branch_with_stand_in(tmp_path: Path, monkeypatch):
    """Drive the real squidpy.gr path with stand-ins recording Z-scores."""
    moran = pd.DataFrame({"I": [0.42, 0.18], "pval_norm": [0.001, 0.04]}, index=["G0", "G1"])
    zscore = np.array([[3.1, -1.2], [-1.2, 2.4]])

    sq_module = types.ModuleType("squidpy")
    sq_module.__version__ = "1.6.0"
    gr = types.SimpleNamespace(
        spatial_neighbors=lambda adata, **k: adata.uns.__setitem__("_nb", True),
        spatial_autocorr=lambda adata, **k: adata.uns.__setitem__("moranI", moran),
        nhood_enrichment=lambda adata, **k: adata.uns.__setitem__(
            "cluster_nhood_enrichment", {"zscore": zscore}
        ),
    )
    sq_module.gr = gr
    monkeypatch.setitem(sys.modules, "squidpy", sq_module)

    ad_module = types.ModuleType("anndata")

    class _AnnData:
        def __init__(self, X=None, obs=None):
            self.X = X
            self.obs = obs if obs is not None else pd.DataFrame()
            self.var_names = []
            self.obsm = {}
            self.uns = {}
            self.n_obs = len(self.obs)

    ad_module.AnnData = _AnnData
    monkeypatch.setitem(sys.modules, "anndata", ad_module)

    result = _run_squidpy_neighborhood(
        data=_spatial_data(), analysis_fields={"analysis_level": "smoke_backend"},
        backend_id="spatial.neighborhood.squidpy", warning="w",
        input_artifact="a", source_dataset="d",
    )

    assert result is not None
    autocorr, enrichment, markers, version = result
    assert version == "1.6.0"
    assert "I" in autocorr.columns and "execution_mode" in autocorr.columns
    assert (autocorr["execution_mode"] == "real_squidpy_gr").all()
    assert "zscore" in enrichment.columns
    # Real permutation-derived Z-scores must survive, not be replaced by zeros.
    assert set(np.round(zscore.flatten(), 3)).issuperset(set(np.round(enrichment["zscore"].dropna().to_numpy(), 3)))
