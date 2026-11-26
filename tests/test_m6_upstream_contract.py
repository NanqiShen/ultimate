"""M6 P0: upstream handoff contract and strict acceptance checker.

The M5 suite could exit 0 with every stage blocked or skipped. These tests pin the
inverted contract: a required check that is missing, skipped, blocked or failed
must always make the run fail, and the generated import config must point at the
paths the Slurm job actually writes.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from ultimate.cli import main
from ultimate.prepare_upstream import (
    ARTIFACTS_SUBDIR,
    prepare_upstream_package,
)
from ultimate.upstream_acceptance import (
    REQUIRED_M6_ROUTES,
    STATUS_BLOCKED,
    STATUS_FAIL,
    STATUS_MISSING,
    STATUS_PASS,
    STATUS_SKIPPED,
    RouteAcceptance,
    check_10x_directory,
    check_downstream_consumed,
    check_expected_files,
    check_matrix_contract,
    check_salmon_quant_reconciliation,
    check_sample_count_agreement,
    evaluate,
    exit_code,
    sha256_file,
    write_acceptance,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# fixtures & test helpers
# --------------------------------------------------------------------------


def _fastq(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x1f\x8b\x08\x00fake")
    return path


def _salmon_index(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "info.json").write_text("{}", encoding="utf-8")
    return root


def _star_index(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "SAindex").write_text("x", encoding="utf-8")
    (root / "Genome").write_text("x", encoding="utf-8")
    return root


def _sheet(tmp_path: Path, header: str = "sample\tfastq_1\tfastq_2", rows: list[str] | None = None) -> Path:
    _fastq(tmp_path / "fastq" / "S1_R1.fastq.gz")
    _fastq(tmp_path / "fastq" / "S1_R2.fastq.gz")
    body = rows if rows is not None else ["S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz"]
    path = tmp_path / "samplesheet.tsv"
    path.write_text("\n".join([header, *body]) + "\n", encoding="utf-8")
    return path


def _package(tmp_path: Path, route: str = "salmon", **kwargs) -> Path:
    index = _salmon_index(tmp_path / "sidx") if route == "salmon" else _star_index(tmp_path / "stidx")
    out = tmp_path / f"pkg_{route}"
    # Never evaluate the default sample sheet eagerly: doing so would rewrite the
    # caller's samplesheet before it is used.
    samplesheet = kwargs.pop("samplesheet", None) or _sheet(tmp_path)
    result = prepare_upstream_package(
        route=route, samplesheet=samplesheet, output_dir=out, index=index, **kwargs,
    )
    assert result.status == "ready", result.blocked_reason
    return out


def _setup_valid_salmon(pkg: Path, sample: str = "S1", *, counts_val: int = 10, num_reads_val: float = 10.2) -> tuple[Path, Path]:
    artifacts = pkg / ARTIFACTS_SUBDIR
    artifacts.mkdir(parents=True, exist_ok=True)
    matrix = artifacts / "counts_raw.tsv"
    matrix.write_text(f"feature_id\t{sample}\nTX1\t{counts_val}\nTX2\t0\n", encoding="utf-8")
    (artifacts / "import_config.yaml").write_text("module: rnaseq\n", encoding="utf-8")

    quant_dir = artifacts / "quant" / sample
    quant_dir.mkdir(parents=True, exist_ok=True)
    (quant_dir / "quant.sf").write_text(
        f"Name\tLength\tEffectiveLength\tTPM\tNumReads\nTX1\t100\t80\t1000\t{num_reads_val}\nTX2\t200\t180\t0\t0.0\n",
        encoding="utf-8",
    )
    mat_sha = sha256_file(matrix)
    manifest = pkg.parent / f"downstream_manifest_{pkg.name}.json"
    manifest.write_text(
        json.dumps({
            "run_id": f"run_rnaseq_{sample}",
            "generated_at": "2026-11-26T12:00:00Z",
            "status": "ready",
            "input_matrix": str(matrix),
            "input_checksum": mat_sha,
            "n_samples": 1,
            "n_genes": 2,
            "delivery_allowed": False,
        }),
        encoding="utf-8",
    )
    return matrix, manifest


def _setup_valid_starsolo(pkg: Path, samples: tuple[str, ...] = ("S1", "S2")) -> tuple[list[Path], Path]:
    artifacts = pkg / ARTIFACTS_SUBDIR
    sample_dirs = []
    matrix_shas = {}
    for s in samples:
        filtered = artifacts / "star" / s / "filtered"
        filtered.mkdir(parents=True, exist_ok=True)
        mtx_content = "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 5\n2 2 3\n"
        with gzip.open(filtered / "matrix.mtx.gz", "wt", encoding="utf-8") as h:
            h.write(mtx_content)
        with gzip.open(filtered / "barcodes.tsv.gz", "wt", encoding="utf-8") as h:
            h.write("AAAC-1\nAAAG-1\n")
        with gzip.open(filtered / "features.tsv.gz", "wt", encoding="utf-8") as h:
            h.write("G1\tG1\tGene Expression\nG2\tG2\tGene Expression\n")
        sample_dirs.append(filtered)
        matrix_shas[s] = sha256_file(filtered / "matrix.mtx.gz")

    (artifacts / "import_config.yaml").write_text("module: scrna\n", encoding="utf-8")
    manifest = pkg.parent / f"downstream_manifest_{pkg.name}.json"
    manifest.write_text(
        json.dumps({
            "run_id": "run_scrna_multi",
            "generated_at": "2026-11-26T12:00:00Z",
            "status": "ready",
            "input_path": str(sample_dirs[0]),
            "samples": [{"sample_id": s, "input_path": str(sd)} for s, sd in zip(samples, sample_dirs)],
            "input_checksums": matrix_shas,
            "n_cells": 4,
            "n_genes": 2,
            "delivery_allowed": False,
        }),
        encoding="utf-8",
    )
    return sample_dirs, manifest


# --------------------------------------------------------------------------
# P0-A: artefact path / import config consistency
# --------------------------------------------------------------------------


def test_salmon_import_config_points_at_actual_artifacts_dir(tmp_path: Path):
    pkg = _package(tmp_path, "salmon")
    config = yaml.safe_load((pkg / "expected_import_project.yaml").read_text(encoding="utf-8"))
    matrix = config["modules"]["rnaseq"]["input_matrix"]
    assert matrix == str(pkg / ARTIFACTS_SUBDIR / "counts_raw.tsv")
    script = (pkg / "run_upstream.sbatch").read_text(encoding="utf-8")
    assert f"OUT_DIR={pkg / ARTIFACTS_SUBDIR}" in script
    assert "counts_raw.tsv" in script


def test_starsolo_import_config_points_at_filtered_per_sample(tmp_path: Path):
    _fastq(tmp_path / "fastq" / "S2_R1.fastq.gz")
    _fastq(tmp_path / "fastq" / "S2_R2.fastq.gz")
    sheet = _sheet(
        tmp_path,
        header="sample\tfastq_1\tfastq_2\tstrandedness",
        rows=["S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz\tforward",
              "S2\tfastq/S2_R1.fastq.gz\tfastq/S2_R2.fastq.gz\tforward"],
    )
    pkg = _package(tmp_path, "starsolo", samplesheet=sheet)
    config = yaml.safe_load((pkg / "expected_import_project.yaml").read_text(encoding="utf-8"))
    samples = config["modules"]["scrna"]["samples"]
    assert [row["sample_id"] for row in samples] == ["S1", "S2"]
    for row in samples:
        assert row["input_path"].endswith(f"star/{row['sample_id']}/filtered")
        assert row["input_path"].startswith(str(pkg / ARTIFACTS_SUBDIR))
    script = (pkg / "run_upstream.sbatch").read_text(encoding="utf-8")
    assert 'filtered="$sample_out/filtered"' in script


def test_count_semantics_are_recorded_honestly(tmp_path: Path):
    salmon = json.loads((_package(tmp_path, "salmon") / "job_package_manifest.json").read_text(encoding="utf-8"))
    sem = salmon["count_semantics"]
    assert sem["level"] == "transcript"
    assert sem["is_gene_level_raw_integer_counts"] == "false"
    assert "NumReads" in sem["unit"]

    stidx_sheet = _sheet(tmp_path / "s2")
    star = json.loads((_package(tmp_path / "s2", "starsolo", samplesheet=stidx_sheet) / "job_package_manifest.json").read_text(encoding="utf-8"))
    assert star["count_semantics"]["level"] == "cell_x_gene"
    assert star["count_semantics"]["is_gene_level_raw_integer_counts"] == "true"


# --------------------------------------------------------------------------
# P0-B: the validator and the job consume one normalised samplesheet
# --------------------------------------------------------------------------


def test_normalized_samplesheet_is_written_and_consumed(tmp_path: Path):
    sheet = _sheet(tmp_path, header="Sample_ID\tFQ1\tFQ2\tStrandedness",
                   rows=["S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz\tforward"])
    pkg = _package(tmp_path, "salmon", samplesheet=sheet)
    normalized = pkg / "samplesheet_normalized.tsv"
    assert normalized.is_file()
    lines = normalized.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "sample_id\tfastq_1\tfastq_2\tstrandedness"
    fields = lines[1].split("\t")
    assert fields[0] == "S1"
    assert fields[1] == "fastq/S1_R1.fastq.gz"
    assert fields[2] == "fastq/S1_R2.fastq.gz"
    assert fields[3] in {"forward", "auto"}
    script = (pkg / "run_upstream.sbatch").read_text(encoding="utf-8")
    assert f"SAMPLESHEET={normalized}" in script
    assert str(sheet) not in script


def test_shuffled_columns_and_aliases_are_normalized(tmp_path: Path):
    sheet = _sheet(tmp_path, header="FQ2\tsample\tFQ1",
                   rows=["fastq/S1_R2.fastq.gz\tS1\tfastq/S1_R1.fastq.gz"])
    pkg = _package(tmp_path, "salmon", samplesheet=sheet)
    row = (pkg / "samplesheet_normalized.tsv").read_text(encoding="utf-8").splitlines()[1].split("\t")
    assert row[0] == "S1"
    assert row[1] == "fastq/S1_R1.fastq.gz"
    assert row[2] == "fastq/S1_R2.fastq.gz"


def test_duplicate_sample_id_blocks(tmp_path: Path):
    sheet = _sheet(tmp_path, rows=[
        "S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz",
        "S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz",
    ])
    result = prepare_upstream_package(
        route="salmon", samplesheet=sheet, output_dir=tmp_path / "pkg",
        index=_salmon_index(tmp_path / "sidx"),
    )
    assert result.status == "blocked"
    assert result.blocked_reason.startswith("duplicate_sample_id")
    assert not (tmp_path / "pkg" / "run_upstream.sbatch").exists()


def test_missing_paired_read_blocks(tmp_path: Path):
    sheet = _sheet(tmp_path, rows=["S1\t-\tfastq/S1_R2.fastq.gz"])
    result = prepare_upstream_package(
        route="salmon", samplesheet=sheet, output_dir=tmp_path / "pkg",
        index=_salmon_index(tmp_path / "sidx"),
    )
    assert result.status == "blocked"
    assert any("missing_paired_read" in blocker for blocker in result.blockers)


def test_paths_with_spaces_are_supported(tmp_path: Path):
    spaced = tmp_path / "fastq with spaces"
    spaced.mkdir()
    _fastq(spaced / "S1 R1.fastq.gz")
    sheet = tmp_path / "samplesheet.tsv"
    sheet.write_text("sample\tfastq_1\nS1\tfastq with spaces/S1 R1.fastq.gz\n", encoding="utf-8")
    result = prepare_upstream_package(
        route="salmon", samplesheet=sheet, output_dir=tmp_path / "pkg with spaces",
        index=_salmon_index(tmp_path / "sidx"),
    )
    assert result.status == "ready", result.blocked_reason
    script = Path(result.artifacts["run_upstream_sbatch"]).read_text(encoding="utf-8")
    assert '"' in script.split("SAMPLESHEET=")[1].splitlines()[0]


def test_existing_non_empty_output_dir_is_refused(tmp_path: Path):
    out = tmp_path / "pkg"
    out.mkdir()
    (out / "previous_result.tsv").write_text("keep me", encoding="utf-8")

    result = prepare_upstream_package(
        route="salmon", samplesheet=_sheet(tmp_path), output_dir=out,
        index=_salmon_index(tmp_path / "sidx"),
    )
    assert result.status == "blocked"
    assert result.blocked_reason.startswith("output_dir_not_empty")
    assert (out / "previous_result.tsv").read_text(encoding="utf-8") == "keep me"
    assert not (out / "run_upstream.sbatch").exists()


# --------------------------------------------------------------------------
# P0-C & P0-D: strict acceptance checker unit contracts
# --------------------------------------------------------------------------


def test_any_non_pass_required_status_fails_the_run():
    for bad in (STATUS_SKIPPED, STATUS_BLOCKED, STATUS_FAIL, STATUS_MISSING):
        route = RouteAcceptance(route="salmon")
        route.add("job_completed", STATUS_PASS)
        route.add("counts_produced", bad, "simulated")
        payload = evaluate([route], required_routes=())
        assert payload["status"] == STATUS_FAIL, bad
        assert exit_code(payload) == 1


def test_all_pass_yields_zero_exit():
    route = RouteAcceptance(route="salmon")
    route.add("job_completed", STATUS_PASS)
    payload = evaluate([route], required_routes=())
    assert payload["status"] == STATUS_PASS
    assert exit_code(payload) == 0


def test_optional_check_failure_does_not_fail_the_run():
    route = RouteAcceptance(route="salmon")
    route.add("job_completed", STATUS_PASS)
    route.add("nice_to_have", STATUS_SKIPPED, required=False)
    assert evaluate([route], required_routes=())["status"] == STATUS_PASS


def test_absent_required_check_is_missing_not_pass():
    route = RouteAcceptance(route="salmon")
    route.add("job_completed", STATUS_PASS)
    route.add("quant_sf_nonempty", STATUS_MISSING, "no quant.sf")
    assert route.failures
    assert evaluate([route], required_routes=())["status"] == STATUS_FAIL


def test_matrix_contract_rejects_bad_domains(tmp_path: Path):
    good = tmp_path / "good.tsv"
    good.write_text("feature_id\tS1\tS2\nG1\t5\t6\nG2\t0\t9\n", encoding="utf-8")
    route = RouteAcceptance(route="salmon")
    check_matrix_contract(route, "m", good, expected_samples=["S1", "S2"], require_integer=True)
    assert not route.failures

    bad = tmp_path / "bad.tsv"
    bad.write_text("feature_id\tS1\tS2\nG1\t1.5\t-2\nG2\t0\t0\n", encoding="utf-8")
    route2 = RouteAcceptance(route="salmon")
    check_matrix_contract(route2, "m", bad, expected_samples=["S1", "S2"], require_integer=True, require_nonzero_signal=True)
    failed = {check.name for check in route2.failures}
    assert any("values_integer" in name for name in failed)
    assert any("values_nonnegative" in name for name in failed)


def test_matrix_sample_mismatch_fails(tmp_path: Path):
    matrix = tmp_path / "m.tsv"
    matrix.write_text("feature_id\tS1\tWRONG\nG1\t5\t6\n", encoding="utf-8")
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample_id\tfastq_1\tfastq_2\tstrandedness\nS1\ta\tb\tforward\n", encoding="utf-8")
    route = RouteAcceptance(route="salmon")
    check_sample_count_agreement(route, "m", matrix, sheet)
    assert route.failures


def test_10x_directory_checks_dims_and_entries(tmp_path: Path):
    good = tmp_path / "filtered"
    good.mkdir()
    (good / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 5\n2 2 6\n", encoding="utf-8")
    (good / "barcodes.tsv").write_text("AAAC-1\nAAAG-1\n", encoding="utf-8")
    (good / "features.tsv").write_text("G1\tG1\tGene Expression\nG2\tG2\tGene Expression\n", encoding="utf-8")
    route = RouteAcceptance(route="starsolo")
    check_10x_directory(route, "s1", good)
    assert not route.failures, [c.name for c in route.failures]

    bad = tmp_path / "bad10x"
    bad.mkdir()
    (bad / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n3 2 1\n1 1 5\n", encoding="utf-8")
    (bad / "barcodes.tsv").write_text("AAAC-1\nAAAG-1\n", encoding="utf-8")
    (bad / "features.tsv").write_text("G1\tG1\tGene Expression\nG2\tG2\tGene Expression\n", encoding="utf-8")
    route2 = RouteAcceptance(route="starsolo")
    check_10x_directory(route2, "s1", bad)
    assert any("dims_consistent" in check.name for check in route2.failures)


# --------------------------------------------------------------------------
# Hermes P0 Regression Tests: 4 Counterexamples and their Valid Controls
# --------------------------------------------------------------------------


def test_m6_counterexample_1_missing_required_route(tmp_path: Path):
    """Counterexample 1: Single route cannot impersonate complete M6 sign-off."""
    pkg_salmon = _package(tmp_path / "s", "salmon")
    _, m_salmon = _setup_valid_salmon(pkg_salmon)

    # 1. Defect case: only salmon provided under default require-both-routes
    runner = CliRunner()
    res_bad = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_salmon), "--downstream-manifest", str(m_salmon)],
    )
    assert res_bad.exit_code != 0
    payload_bad = json.loads(res_bad.output)
    assert payload_bad["status"] == "fail"
    assert "starsolo" in payload_bad["missing_routes"]

    # 2. Valid control: both salmon and starsolo routes provided
    pkg_star = _package(tmp_path / "st", "starsolo")
    _, m_star = _setup_valid_starsolo(pkg_star, samples=("S1",))
    res_good = runner.invoke(
        main,
        [
            "check-upstream-acceptance",
            "--package-dir", str(pkg_salmon), "--downstream-manifest", str(m_salmon),
            "--package-dir", str(pkg_star), "--downstream-manifest", str(m_star),
        ],
    )
    assert res_good.exit_code == 0, res_good.output
    payload_good = json.loads(res_good.output)
    assert payload_good["status"] == "pass"
    assert payload_good["missing_routes"] == []


def test_m6_counterexample_2_salmon_quant_reconciliation(tmp_path: Path):
    """Counterexample 2: Mismatched, missing, or non-finite quant.sf fails."""
    runner = CliRunner()

    # 2a: Mismatch between quant.sf and counts_raw.tsv
    pkg_mismatch = _package(tmp_path / "mismatch", "salmon")
    _, m_path = _setup_valid_salmon(pkg_mismatch, counts_val=99, num_reads_val=10.0)
    res_mismatch = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_mismatch), "--downstream-manifest", str(m_path), "--allow-single-route"],
    )
    assert res_mismatch.exit_code != 0
    assert "quant_reconciliation" in res_mismatch.output

    # 2b: Missing quant.sf
    pkg_missing = _package(tmp_path / "missing", "salmon")
    _, m_miss_path = _setup_valid_salmon(pkg_missing)
    (pkg_missing / ARTIFACTS_SUBDIR / "quant" / "S1" / "quant.sf").unlink()
    res_missing = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_missing), "--downstream-manifest", str(m_miss_path), "--allow-single-route"],
    )
    assert res_missing.exit_code != 0
    assert "quant_sf_present" in res_missing.output

    # 2c: Non-finite NaN/Inf in matrix
    pkg_nan = _package(tmp_path / "nan", "salmon")
    mat_nan, m_nan_path = _setup_valid_salmon(pkg_nan)
    mat_nan.write_text("feature_id\tS1\nTX1\tNaN\nTX2\t0\n", encoding="utf-8")
    res_nan = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_nan), "--downstream-manifest", str(m_nan_path), "--allow-single-route"],
    )
    assert res_nan.exit_code != 0
    assert "values_finite" in res_nan.output or "quant_reconciliation" in res_nan.output

    # Valid control: exact reconciliation and finite values
    pkg_ctrl = _package(tmp_path / "ctrl", "salmon")
    _, m_ctrl_path = _setup_valid_salmon(pkg_ctrl, counts_val=10, num_reads_val=10.0)
    res_ctrl = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_ctrl), "--downstream-manifest", str(m_ctrl_path), "--allow-single-route"],
    )
    assert res_ctrl.exit_code == 0, res_ctrl.output
    assert json.loads(res_ctrl.output)["status"] == "pass"


def test_m6_counterexample_3_starsolo_matrix_entries(tmp_path: Path):
    """Counterexample 3: Matrix Market entry corruptions or missing samples fail."""
    runner = CliRunner()

    # 3a: nnz declared 5 but only 1 data line
    pkg_nnz = _package(tmp_path / "nnz", "starsolo")
    sample_dirs, m_nnz = _setup_valid_starsolo(pkg_nnz, samples=("S1",))
    with gzip.open(sample_dirs[0] / "matrix.mtx.gz", "wt", encoding="utf-8") as h:
        h.write("%%MatrixMarket matrix coordinate integer general\n2 2 5\n1 1 5\n")
    res_nnz = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_nnz), "--downstream-manifest", str(m_nnz), "--allow-single-route"],
    )
    assert res_nnz.exit_code != 0
    assert "entry_count_matches_nnz" in res_nnz.output

    # 3b: Entry has negative or float count
    pkg_badval = _package(tmp_path / "badval", "starsolo")
    sdirs_badval, m_badval = _setup_valid_starsolo(pkg_badval, samples=("S1",))
    with gzip.open(sdirs_badval[0] / "matrix.mtx.gz", "wt", encoding="utf-8") as h:
        h.write("%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 -5\n2 2 3.14\n")
    res_badval = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_badval), "--downstream-manifest", str(m_badval), "--allow-single-route"],
    )
    assert res_badval.exit_code != 0
    assert "values_integer_nonnegative" in res_badval.output

    # 3c: Nonzero signal missing (all counts zero)
    pkg_zero = _package(tmp_path / "zero", "starsolo")
    sdirs_zero, m_zero = _setup_valid_starsolo(pkg_zero, samples=("S1",))
    with gzip.open(sdirs_zero[0] / "matrix.mtx.gz", "wt", encoding="utf-8") as h:
        h.write("%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 0\n2 2 0\n")
    res_zero = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_zero), "--downstream-manifest", str(m_zero), "--allow-single-route"],
    )
    assert res_zero.exit_code != 0
    assert "nonzero_signal" in res_zero.output

    # 3d: Multi-sample missing S2 filtered dir
    sheet_multi = _sheet(
        tmp_path / "multi",
        header="sample\tfastq_1\tfastq_2",
        rows=["S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz",
              "S2\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz"],
    )
    pkg_multi = _package(tmp_path / "multi", "starsolo", samplesheet=sheet_multi)
    # Only create S1
    sdirs_multi, m_multi = _setup_valid_starsolo(pkg_multi, samples=("S1",))
    res_multi = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_multi), "--downstream-manifest", str(m_multi), "--allow-single-route"],
    )
    assert res_multi.exit_code != 0
    assert "upstream:S2" in res_multi.output

    # Valid control: multi-sample all valid
    sheet_good = _sheet(
        tmp_path / "good",
        header="sample\tfastq_1\tfastq_2",
        rows=["S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz",
              "S2\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz"],
    )
    pkg_good = _package(tmp_path / "good", "starsolo", samplesheet=sheet_good)
    _setup_valid_starsolo(pkg_good, samples=("S1", "S2"))
    m_good_path = tmp_path / "good" / f"downstream_manifest_{pkg_good.name}.json"
    res_good = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_good), "--downstream-manifest", str(m_good_path), "--allow-single-route"],
    )
    assert res_good.exit_code == 0, res_good.output
    assert json.loads(res_good.output)["status"] == "pass"


def test_m6_counterexample_4_downstream_manifest_evidence(tmp_path: Path):
    """Counterexample 4: Downstream must provide runtime binding, checksum, dims, and delivery_allowed=False."""
    runner = CliRunner()
    pkg = _package(tmp_path, "salmon")
    mat, m_valid = _setup_valid_salmon(pkg)
    valid_payload = json.loads(m_valid.read_text(encoding="utf-8"))

    # 4a: Missing runtime binding (missing run_id or generated_at)
    p_no_runid = dict(valid_payload)
    del p_no_runid["run_id"]
    m_no_runid = tmp_path / "m_no_runid.json"
    m_no_runid.write_text(json.dumps(p_no_runid), encoding="utf-8")
    res_4a = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_no_runid), "--allow-single-route"],
    )
    assert res_4a.exit_code != 0
    assert "downstream_runtime_binding" in res_4a.output

    # 4b: Checksum missing or mismatched
    p_bad_sha = dict(valid_payload, input_checksum="0123456789abcdef" * 4)
    m_bad_sha = tmp_path / "m_bad_sha.json"
    m_bad_sha.write_text(json.dumps(p_bad_sha), encoding="utf-8")
    res_4b = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_bad_sha), "--allow-single-route"],
    )
    assert res_4b.exit_code != 0
    assert "downstream_input_checksum" in res_4b.output

    # 4c: Dimensions missing
    p_no_dims = dict(valid_payload)
    del p_no_dims["n_samples"]
    del p_no_dims["n_genes"]
    m_no_dims = tmp_path / "m_no_dims.json"
    m_no_dims.write_text(json.dumps(p_no_dims), encoding="utf-8")
    res_4c = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_no_dims), "--allow-single-route"],
    )
    assert res_4c.exit_code != 0
    assert "downstream_actual_dimensions" in res_4c.output

    # 4d: delivery_allowed missing, string "false", or True
    for bad_delivery in (None, "false", True):
        p_bad_gate = dict(valid_payload)
        if bad_delivery is None:
            del p_bad_gate["delivery_allowed"]
        else:
            p_bad_gate["delivery_allowed"] = bad_delivery
        m_bad_gate = tmp_path / f"m_bad_gate_{bad_delivery}.json"
        m_bad_gate.write_text(json.dumps(p_bad_gate), encoding="utf-8")
        res_4d = runner.invoke(
            main,
            ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_bad_gate), "--allow-single-route"],
        )
        assert res_4d.exit_code != 0, f"failed to reject bad delivery_allowed={bad_delivery}"
        assert "delivery_allowed_boolean_false" in res_4d.output

    # Valid control: valid runtime binding, matched checksum, dims, and delivery_allowed=False
    res_ctrl = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_valid), "--allow-single-route"],
    )
    assert res_ctrl.exit_code == 0, res_ctrl.output
    assert json.loads(res_ctrl.output)["status"] == "pass"


# --------------------------------------------------------------------------
# templates & auxiliary suite contract
# --------------------------------------------------------------------------


def test_templates_require_normalised_samplesheet():
    for name in ("slurm/upstream_bulk_rnaseq_salmon.sbatch", "slurm/upstream_scrna_starsolo.sbatch"):
        text = (REPO_ROOT / name).read_text(encoding="utf-8")
        assert "samplesheet_normalized" in text or "normalised table" in text, name
        assert "sample_id\\tfastq_1\\tfastq_2\\tstrandedness" in text, name


def test_auxiliary_suite_fails_on_non_ready_stage():
    text = (REPO_ROOT / "slurm" / "xielab_e2e_regression_suite.sbatch").read_text(encoding="utf-8")
    # The old too-loose check must be gone.
    assert 'grep -q $\'\\tfailed\\t\'' not in text
    assert 'not_ready=' in text
    assert "exit 1" in text
    assert "ULTIMATE_E2E_ALLOW_SKIP" in text
