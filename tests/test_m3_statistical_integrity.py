"""M3 scientific-rigor regression tests.

Pins the four claims of Phase 3:

* KM/log-rank come from real data (lifelines/R), not a fixed exponential curve,
  and a shuffled clinical table still gives the correct paired result;
* Welch t does not inflate significance on small, unequal-variance cohorts;
* the HTML report embeds real <img> grids and Top-10 tables, and the customer
  package keeps image paths valid while redacting internal paths;
* no Python backend writes JSON text into a .rds/.RData filename.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ultimate import survival_backend
from ultimate.bulk import _two_group_stats, _wgcna_tables
from ultimate.customer_package import build_customer_package
from ultimate.demo import init_project
from ultimate.pipeline import run_pipeline_from_config
from ultimate.report import build_report

REPO_ROOT = Path(__file__).resolve().parents[1]


# --- 1. survival analysis ----------------------------------------------------


def test_no_fake_exponential_curves_in_source():
    source_files = list((REPO_ROOT / "src").rglob("*.py"))
    offenders = [str(path) for path in source_files if "np.exp(-months" in path.read_text(encoding="utf-8")]
    assert offenders == []


def test_shuffled_samples_still_produce_correct_km_and_logrank(tmp_path: Path):
    """A shuffled clinical table must not break sample alignment."""
    n = 30
    ids = [f"S{i}" for i in range(n)]
    # High group has clearly shorter survival; log-rank should be significant.
    time = np.array([6.0] * 15 + [40.0] * 15)
    event = np.array([1] * 15 + [0] * 15)
    group = np.array(["high"] * 15 + ["low"] * 15)
    order = np.random.default_rng(42).permutation(n)
    clinical = pd.DataFrame({"sample_id": [ids[i] for i in order],
                             "survival_time": time[order], "event": event[order]})
    expression = pd.Series([2.0 if g == "high" else 0.0 for g in group], index=ids)
    result = survival_backend.run_survival_analysis(
        clinical, expression,
        figures_dir=tmp_path / "fig", tables_dir=tmp_path / "tab", label="GENE",
    )
    assert result.status == "ready", result.reason
    assert result.engine in {"lifelines", "r_survival"}
    assert result.logrank_p is not None and result.logrank_p < 0.05
    assert result.n_samples == n
    assert Path(result.figure_path).is_file()
    assert not result.km_table.empty


def test_missing_survival_columns_returns_explicit_skipped(tmp_path: Path):
    clinical = pd.DataFrame({"sample_id": ["S1", "S2"], "age": [50, 60]})
    result = survival_backend.run_survival_analysis(
        clinical, pd.Series([1.0, 0.0], index=["S1", "S2"]),
        figures_dir=tmp_path / "fig", tables_dir=tmp_path / "tab", label="GENE",
    )
    assert result.status == "skipped"
    assert "missing_survival_columns" in result.reason
    assert Path(result.figure_path).is_file()


def test_survival_engine_absence_is_explicit_not_fabricated(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(survival_backend, "_available_engine", lambda module: "")
    clinical = pd.DataFrame({"sample_id": ["S1", "S2", "S3", "S4"],
                             "survival_time": [1.0, 2.0, 3.0, 4.0], "event": [1, 1, 0, 1]})
    result = survival_backend.run_survival_analysis(
        clinical, pd.Series([2, 2, 0, 0], index=clinical["sample_id"]),
        figures_dir=tmp_path / "fig", tables_dir=tmp_path / "tab", label="G",
    )
    assert result.status == "skipped"
    assert "no_survival_engine" in result.reason


# --- 2. Welch t --------------------------------------------------------------


def test_welch_t_small_sample_is_not_inflated():
    """Unequal variances with n=4 per group: Welch p must be >= naive z p."""
    rng = np.random.default_rng(11)
    ctrl = rng.normal(5, 0.5, (20, 4))
    case = rng.normal(5.6, 4.0, (20, 4))  # same mean shift, much larger variance
    matrix = pd.DataFrame(
        np.concatenate([ctrl, case], axis=1),
        index=[f"G{i}" for i in range(20)],
        columns=[f"C{i}" for i in range(4)] + [f"T{i}" for i in range(4)],
    )
    result = _two_group_stats(matrix, list(matrix.columns[:4]), list(matrix.columns[4:]))
    assert (result["test_method"] == "welch_t_scipy").all()
    assert "df" in result.columns
    # Welch df is well below the pooled df for unequal variances.
    assert result["df"].max() < 7
    # At least some features must not be declared significant despite the shift.
    assert (result["pvalue"] > 0.05).any()


def test_welch_bh_fdr_is_monotonic():
    rng = np.random.default_rng(5)
    matrix = pd.DataFrame(rng.normal(size=(40, 8)), index=[f"G{i}" for i in range(40)])
    result = _two_group_stats(matrix, list(matrix.columns[:4]), list(matrix.columns[4:]))
    assert (result["padj"].diff().dropna() >= -1e-9).all()
    assert result["pvalue"].between(0, 1).all()


# --- WGCNA is never a quartile proxy ----------------------------------------


def test_wgcna_reports_skip_not_four_quantile_bins(tmp_path: Path):
    matrix = pd.DataFrame(np.random.default_rng(1).normal(size=(30, 10)),
                          index=[f"G{i}" for i in range(30)], columns=[f"S{i}" for i in range(10)])
    clinical = pd.DataFrame({"sample_id": list(matrix.columns)})
    tables = tmp_path / "tables"; tables.mkdir()
    result = _wgcna_tables(matrix, clinical, tables)
    assert result["wgcna_status"] == "skipped"
    assert "insufficient_samples" in result["wgcna_skip_reason"]
    frame = pd.read_csv(tables / "wgcna_module_assignments.tsv", sep="\t")
    assert "blue" not in set(frame["module"].astype(str))
    assert "skip_reason" in frame.columns


# --- 3. HTML report ----------------------------------------------------------


def _synth_run(tmp_path: Path) -> Path:
    project = init_project("rnaseq", tmp_path / "proj", demo_data=True)
    manifest = run_pipeline_from_config(Path(project["config_path"]))
    run_dir = Path(manifest["output_dir"])
    build_report(run_dir)
    return run_dir


def test_report_html_contains_image_grid_and_top_table(tmp_path: Path):
    run_dir = _synth_run(tmp_path)
    html = (run_dir / "reports" / "report.html").read_text(encoding="utf-8")
    assert "<img" in html
    assert "figure-grid" in html
    assert "Top 10" in html
    assert 'class="card"' in html
    assert "解释边界" in html


def test_report_embed_images_uses_base64(tmp_path: Path):
    run_dir = _synth_run(tmp_path)
    manifest_path = run_dir / "run_manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["report"] = {"embed_images": True}
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    build_report(run_dir)
    html = (run_dir / "reports" / "report.html").read_text(encoding="utf-8")
    assert 'src="data:image/png;base64,' in html


def test_customer_package_rewrites_image_paths_and_sanitizes(tmp_path: Path):
    run_dir = _synth_run(tmp_path)
    report = run_dir / "reports" / "report.html"
    report.write_text(report.read_text(encoding="utf-8") + "<p>/dshare/xielab/rawdata/TFC/shennq</p>", encoding="utf-8")
    package = tmp_path / "customer"
    build_customer_package(run_dir=run_dir, output_dir=package)
    html = (package / "report.html").read_text(encoding="utf-8")
    srcs = re.findall(r'<img[^>]*src="([^"]+)"', html)
    assert srcs, "customer report must render images"
    for src in srcs:
        assert not src.startswith("/")
        assert (package / src).exists(), f"broken customer image path: {src}"
    assert "/dshare" not in html and "/shared" not in html


# --- 4. no Python-written JSON disguised as .rds/.RData ----------------------


def test_no_python_source_writes_json_into_rds_suffix():
    offenders = []
    for path in (REPO_ROOT / "src" / "ultimate").rglob("*.py"):
        lines = path.read_text(encoding="utf-8").splitlines()
        for number, line in enumerate(lines, 1):
            if "write_text" not in line or "json" not in line:
                continue
            window = "\n".join(lines[max(0, number - 40):number])
            if re.search(r'\.(rds|RData)"', window):
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{number}")
    assert offenders == []


def test_real_r_scripts_still_write_rds_via_saverds():
    """The cleanup must not remove genuine R saveRDS outputs."""
    r_dir = REPO_ROOT / "scripts" / "R"
    assert (r_dir / "rnaseq_de_backend.R").read_text().count("saveRDS") >= 1
    assert (r_dir / "survival_backend.R").exists()


def test_survival_r_script_is_packaged():
    from ultimate.env_router import resolve_packaged_asset

    assert resolve_packaged_asset("scripts", "R", "survival_backend.R").is_file()
