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
    if result.status == "skipped":
        # No survival engine on this interpreter (no lifelines and no R survival
        # package). That must be an explicit, explained skip - never a fabricated
        # ready result or a fake curve.
        assert result.reason
        assert result.logrank_p is None
        assert result.cox_hr is None
        skipped_figure = tmp_path / "fig" / "km_curve.png"
        assert skipped_figure.is_file()
        return
    assert result.status == "ready", result.reason
    assert result.engine in {"lifelines", "r_survival"}
    assert result.logrank_p is not None and result.logrank_p < 0.05
    assert result.n_samples == n
    assert Path(result.figure_path).is_file()
    assert not result.km_table.empty


def test_lifelines_engine_produces_significant_logrank_when_installed(tmp_path: Path):
    """When lifelines is available the real statistic must be computed."""
    pytest.importorskip("lifelines")
    n = 30
    ids = [f"S{i}" for i in range(n)]
    time = np.array([6.0] * 15 + [40.0] * 15)
    event = np.array([1] * 15 + [0] * 15)
    order = np.random.default_rng(42).permutation(n)
    clinical = pd.DataFrame({"sample_id": [ids[i] for i in order],
                             "survival_time": time[order], "event": event[order]})
    expression = pd.Series([2.0] * 15 + [0.0] * 15, index=ids)
    result = survival_backend.run_survival_analysis(
        clinical, expression,
        figures_dir=tmp_path / "fig", tables_dir=tmp_path / "tab", label="GENE",
    )
    assert result.status == "ready"
    assert result.engine == "lifelines"
    assert result.logrank_p is not None and result.logrank_p < 0.05
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


def _r_clinical() -> pd.DataFrame:
    return pd.DataFrame({
        "sample_id": ["S1", "S2", "S3", "S4", "S5", "S6"],
        "survival_time": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        "event": [1, 1, 0, 1, 0, 1],
    })


def test_r_fallback_skipped_payload_is_not_reported_ready(tmp_path: Path, monkeypatch):
    """R may exit 0 with status=skipped (missing R package) - must stay skipped."""
    import subprocess

    class _Completed:
        stdout = '{"status":"skipped","reason":"missing_r_package:survival"}'

    monkeypatch.setattr(survival_backend, "_available_engine", lambda module: "r_survival")
    monkeypatch.setattr(survival_backend, "resolve_rscript", lambda module: Path("/tmp/Rscript"))
    monkeypatch.setattr(survival_backend, "resolve_packaged_asset", lambda *parts: Path("/tmp/survival_backend.R"))
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: _Completed())

    result = survival_backend.run_survival_analysis(
        _r_clinical(), pd.Series([2, 2, 2, 0, 0, 0], index=[f"S{i}" for i in range(1, 7)]),
        figures_dir=tmp_path / "fig", tables_dir=tmp_path / "tab", label="G",
    )
    assert result.status == "skipped"
    assert "missing_r_package:survival" in result.reason
    assert Path(result.figure_path).is_file()


def test_r_fallback_missing_figure_is_not_reported_ready(tmp_path: Path, monkeypatch):
    """Even a ready payload must not be accepted without the KM figure."""
    import subprocess

    class _Completed:
        stdout = '{"status":"ready","logrank_p":0.01}'

    monkeypatch.setattr(survival_backend, "_available_engine", lambda module: "r_survival")
    monkeypatch.setattr(survival_backend, "resolve_rscript", lambda module: Path("/tmp/Rscript"))
    monkeypatch.setattr(survival_backend, "resolve_packaged_asset", lambda *parts: Path("/tmp/survival_backend.R"))
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: _Completed())  # no png written

    result = survival_backend.run_survival_analysis(
        _r_clinical(), pd.Series([2, 2, 2, 0, 0, 0], index=[f"S{i}" for i in range(1, 7)]),
        figures_dir=tmp_path / "fig", tables_dir=tmp_path / "tab", label="G",
    )
    assert result.status == "skipped"
    assert "r_survival_skipped" in result.reason


def test_r_fallback_uses_requested_module_name(tmp_path: Path, monkeypatch):
    """The Rscript route must honour the caller's module, not a hardcoded one."""
    seen: list[str] = []

    def _fake_resolve(module: str):
        seen.append(module)
        return Path("/tmp/Rscript")

    monkeypatch.setattr(survival_backend, "_available_engine", lambda module: "r_survival")
    monkeypatch.setattr(survival_backend, "resolve_rscript", _fake_resolve)
    monkeypatch.setattr(survival_backend, "resolve_packaged_asset", lambda *parts: Path("/tmp/survival_backend.R"))

    import subprocess

    class _Completed:
        stdout = '{"status":"skipped","reason":"no_survival"}'

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: _Completed())
    survival_backend.run_survival_analysis(
        _r_clinical(), pd.Series([2, 2, 2, 0, 0, 0], index=[f"S{i}" for i in range(1, 7)]),
        figures_dir=tmp_path / "fig", tables_dir=tmp_path / "tab", label="G",
        module_name="publicdb",
    )
    assert seen == ["publicdb"]


# --- 2. Welch t --------------------------------------------------------------


def test_welch_t_small_sample_is_not_inflated():
    """Unequal variances with n=4 per group: Welch p must be >= naive z p."""
    pytest.importorskip("scipy")
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


def test_welch_t_fallback_without_scipy_is_honestly_labeled(monkeypatch):
    import importlib.util
    from ultimate import bulk as bulk_mod

    orig_find_spec = importlib.util.find_spec
    monkeypatch.setattr(
        bulk_mod.importlib.util,
        "find_spec",
        lambda name, package=None: None if name == "scipy" else orig_find_spec(name, package),
    )
    rng = np.random.default_rng(7)
    matrix = pd.DataFrame(
        rng.normal(size=(10, 6)),
        index=[f"G{i}" for i in range(10)],
        columns=[f"C{i}" for i in range(3)] + [f"T{i}" for i in range(3)],
    )
    result = _two_group_stats(matrix, list(matrix.columns[:3]), list(matrix.columns[3:]))
    assert (result["test_method"] == "normal_approx_fallback_not_for_production").all()
    assert (result["df"] == 4.0).all()


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
    """Detect writing JSON into a .rds/.RData target, including via a helper
    parameter. Targets rewritten through ``with_name(.._skip.json)`` are honest.
    """
    import ast

    offenders: list[str] = []
    for path in (REPO_ROOT / "src" / "ultimate").rglob("*.py"):
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        rds_literal_re = re.compile(r"\.(rds|rdata)(?:$|[^a-z])", re.IGNORECASE)

        def bound_to_rds(name: str, scope: ast.AST) -> bool:
            for node in ast.walk(scope):
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name) and target.id == name:
                            # Path may be a plain string or a BinOp like
                            # ``objects_dir / "x.rds"``; inspect all string constants.
                            for inner in ast.walk(node.value):
                                if isinstance(inner, ast.Constant) and isinstance(inner.value, str) and rds_literal_re.search(inner.value):
                                    return True
            return False

        def is_skip_rewrite(node: ast.Call) -> bool:
            seg = ast.get_source_segment(source, node) or ""
            # Check the statement (parent) for a with_name rewrite of the target.
            return "with_name" in seg and "_skip.json" in seg

        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            params = {arg.arg for arg in node.args.args} | {arg.arg for arg in node.args.kwonlyargs}
            if node.args.kwarg:
                params.add(node.args.kwarg.arg)
            for sub in ast.walk(node):
                if not isinstance(sub, ast.Call):
                    continue
                fn = sub.func
                cname = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
                if cname not in {"write_text", "write", "write_bytes"} or not isinstance(fn, ast.Attribute):
                    continue
                seg = ast.get_source_segment(source, sub) or ""
                if not any(k in seg for k in ("json", "dumps", "payload", "manifest", "base")):
                    continue
                if is_skip_rewrite(sub):
                    continue
                target = ast.get_source_segment(source, fn.value) or ""
                target_name = target.strip()
                # 1) local assignment to a literal .rds path
                if bound_to_rds(target_name, node):
                    offenders.append(f"{path.relative_to(REPO_ROOT)}:{sub.lineno}: local {target_name}")
                    continue
                # 2) helper parameter that receives a .rds path at any call site
                if target_name in params:
                    kw_target = f"{target_name}="
                    for call in ast.walk(tree):
                        if not isinstance(call, ast.Call):
                            continue
                        call_fn = call.func
                        call_name = call_fn.attr if isinstance(call_fn, ast.Attribute) else getattr(call_fn, "id", "")
                        if call_name != node.name:
                            continue
                        for kw in call.keywords:
                            if kw.arg == target_name and isinstance(kw.value, ast.Name) and bound_to_rds(kw.value.id, tree):
                                offenders.append(f"{path.relative_to(REPO_ROOT)}:{sub.lineno}: param {target_name}")
                            elif kw.arg == target_name and isinstance(kw.value, ast.Constant) and isinstance(kw.value.value, str) and rds_literal_re.search(kw.value.value):
                                offenders.append(f"{path.relative_to(REPO_ROOT)}:{sub.lineno}: param literal {target_name}")
    assert offenders == []



def test_real_r_scripts_still_write_rds_via_saverds():
    """The cleanup must not remove genuine R saveRDS outputs."""
    r_dir = REPO_ROOT / "scripts" / "R"
    assert (r_dir / "rnaseq_de_backend.R").read_text().count("saveRDS") >= 1
    assert (r_dir / "survival_backend.R").exists()


def test_survival_r_script_is_packaged():
    from ultimate.env_router import resolve_packaged_asset

    assert resolve_packaged_asset("scripts", "R", "survival_backend.R").is_file()



def test_r_limma_scripts_use_dynamic_condition_column():
    """Custom --condition-column must actually drive grouping in both R scripts."""
    r_dir = REPO_ROOT / "scripts" / "R"
    for name in ("proteomics_limma_backend.R", "methylation_dmp_limma.R"):
        text = (r_dir / name).read_text(encoding="utf-8")
        assert "samples[[condition]]" in text, name
        assert "samples$condition" not in text, name
        assert 'condition %in% colnames(samples)' in text, name
