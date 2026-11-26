"""M6 P0: upstream handoff contract and strict acceptance checker.

The M5 suite could exit 0 with every stage blocked or skipped. These tests pin the
inverted contract: a required check that is missing, skipped, blocked or failed
must always make the run fail, and the generated import config must point at the
paths the Slurm job actually writes.
"""

from __future__ import annotations

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
    STATUS_BLOCKED,
    STATUS_FAIL,
    STATUS_MISSING,
    STATUS_PASS,
    STATUS_SKIPPED,
    RouteAcceptance,
    check_matrix_contract,
    check_sample_count_agreement,
    check_10x_directory,
    evaluate,
    exit_code,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# fixtures
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


# --------------------------------------------------------------------------
# P0-A: artefact path / import config consistency
# --------------------------------------------------------------------------


def test_salmon_import_config_points_at_actual_artifacts_dir(tmp_path: Path):
    pkg = _package(tmp_path, "salmon")
    config = yaml.safe_load((pkg / "expected_import_project.yaml").read_text(encoding="utf-8"))
    matrix = config["modules"]["rnaseq"]["input_matrix"]
    assert matrix == str(pkg / ARTIFACTS_SUBDIR / "counts_raw.tsv")
    # The path must be exactly where the sbatch template writes.
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
        # Must be the filtered/ directory itself, never its parent.
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
    # Strandedness is passed through when supplied, defaulting to "auto".
    assert fields[3] in {"forward", "auto"}
    # The generated job must read the normalised table, not the raw one.
    script = (pkg / "run_upstream.sbatch").read_text(encoding="utf-8")
    assert f"SAMPLESHEET={normalized}" in script
    assert str(sheet) not in script


def test_shuffled_columns_and_aliases_are_normalized(tmp_path: Path):
    """Column order and case must not matter as long as aliases are known."""
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
    # Paths containing spaces must be quoted in the generated script.
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
    # Never overwrite history.
    assert (out / "previous_result.tsv").read_text(encoding="utf-8") == "keep me"
    assert not (out / "run_upstream.sbatch").exists()


# --------------------------------------------------------------------------
# P0-C: strict acceptance checker
# --------------------------------------------------------------------------


def test_any_non_pass_required_status_fails_the_run():
    for bad in (STATUS_SKIPPED, STATUS_BLOCKED, STATUS_FAIL, STATUS_MISSING):
        route = RouteAcceptance(route="salmon")
        route.add("job_completed", STATUS_PASS)
        route.add("counts_produced", bad, "simulated")
        payload = evaluate([route])
        assert payload["status"] == STATUS_FAIL, bad
        assert exit_code(payload) == 1
        assert "counts_produced" in payload["failed_routes"] or payload["failed_routes"] == ["salmon"]


def test_all_pass_yields_zero_exit():
    route = RouteAcceptance(route="salmon")
    route.add("job_completed", STATUS_PASS)
    payload = evaluate([route])
    assert payload["status"] == STATUS_PASS
    assert exit_code(payload) == 0


def test_optional_check_failure_does_not_fail_the_run():
    route = RouteAcceptance(route="salmon")
    route.add("job_completed", STATUS_PASS)
    route.add("nice_to_have", STATUS_SKIPPED, required=False)
    assert evaluate([route])["status"] == STATUS_PASS


def test_absent_required_check_is_missing_not_pass():
    """A check that was never evaluated must not be treated as passing."""
    route = RouteAcceptance(route="salmon")
    route.add("job_completed", STATUS_PASS)
    route.add("quant_sf_nonempty", STATUS_MISSING, "no quant.sf")
    assert route.failures
    assert evaluate([route])["status"] == STATUS_FAIL


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


def test_10x_directory_checks_dims(tmp_path: Path):
    import gzip

    good = tmp_path / "filtered"
    good.mkdir()
    (good / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 5\n2 2 6\n1 2 7\n", encoding="utf-8")
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


def test_cli_acceptance_exits_nonzero_when_upstream_outputs_missing(tmp_path: Path):
    pkg = _package(tmp_path, "salmon")
    manifest = tmp_path / "run_manifest.json"
    manifest.write_text(json.dumps({"status": "ready", "delivery_allowed": False}), encoding="utf-8")

    result = CliRunner().invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(manifest)],
    )

    assert result.exit_code != 0  # counts_raw.tsv was never produced
    payload = json.loads(result.output)
    assert payload["status"] == "fail"


def test_cli_acceptance_passes_with_real_artifacts(tmp_path: Path):
    pkg = _package(tmp_path, "salmon")
    artifacts = pkg / ARTIFACTS_SUBDIR
    artifacts.mkdir(parents=True, exist_ok=True)
    matrix = artifacts / "counts_raw.tsv"
    matrix.write_text("feature_id\tS1\nTX1\t10\nTX2\t0\n", encoding="utf-8")
    (artifacts / "import_config.yaml").write_text("module: rnaseq\n", encoding="utf-8")
    manifest = tmp_path / "run_manifest.json"
    manifest.write_text(
        json.dumps({"status": "ready", "input_matrix": str(matrix), "n_samples": 1, "delivery_allowed": False}),
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(manifest)],
    )

    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["status"] == "pass"


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
