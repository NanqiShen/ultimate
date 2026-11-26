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
from ultimate.pipeline import run_pipeline_from_config
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
            "input_checksums": {matrix.name: mat_sha, str(matrix.resolve()): mat_sha},
            "n_samples": 1,
            "n_features": 2,
            "n_genes": 2,
            "raw_dimensions": {"n_features": 2, "n_samples": 1},
            "delivery_allowed": False,
        }),
        encoding="utf-8",
    )
    return matrix, manifest


def _setup_valid_starsolo(pkg: Path, samples: tuple[str, ...] = ("S1", "S2")) -> tuple[list[Path], Path]:
    artifacts = pkg / ARTIFACTS_SUBDIR
    sample_dirs = []
    all_triad_shas: dict[str, str] = {}
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
        for stem in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"):
            sha = sha256_file(filtered / stem)
            all_triad_shas[stem] = sha
            all_triad_shas[str(filtered / stem)] = sha
            all_triad_shas[f"{s}/{stem}"] = sha

    (artifacts / "import_config.yaml").write_text("module: scrna\n", encoding="utf-8")
    manifest = pkg.parent / f"downstream_manifest_{pkg.name}.json"
    manifest.write_text(
        json.dumps({
            "run_id": "run_scrna_multi",
            "generated_at": "2026-11-26T12:00:00Z",
            "status": "ready",
            "input_path": str(sample_dirs[0]),
            "samples": [{"sample_id": s, "input_path": str(sd)} for s, sd in zip(samples, sample_dirs)],
            "input_checksums": all_triad_shas,
            "n_cells": 2,
            "n_genes": 2,
            "raw_dimensions": {"n_cells": 2, "n_genes": 2},
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
    assert lines[0].startswith("sample_id\tfastq_1\tfastq_2\tstrandedness")
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
# Codex P0 Final Tightening: Real Ingestion Evidence & Four Counterexamples
# --------------------------------------------------------------------------


def test_real_local_ingestion_produces_accepted_manifests(tmp_path: Path):
    """Positive test: Downstream manifests must be produced by REAL reading code."""
    runner = CliRunner()

    # 1. Real Salmon/Bulk ingestion: run_pipeline_from_config
    pkg_salmon = _package(tmp_path / "salmon_real", "salmon")
    art_salmon = pkg_salmon / ARTIFACTS_SUBDIR
    art_salmon.mkdir(parents=True, exist_ok=True)
    matrix_bulk = art_salmon / "counts_raw.tsv"
    matrix_bulk.write_text("feature_id\tS1\nG1\t10\nG2\t20\n", encoding="utf-8")
    (art_salmon / "import_config.yaml").write_text("module: rnaseq\n", encoding="utf-8")
    q_dir = art_salmon / "quant" / "S1"
    q_dir.mkdir(parents=True, exist_ok=True)
    (q_dir / "quant.sf").write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\nG1\t100\t80\t100\t10.0\nG2\t200\t180\t200\t20.0\n",
        encoding="utf-8",
    )
    # Real run: pipeline reads matrix_bulk, calculates sha256, binds run_id and raw_dimensions
    bulk_run_manifest = run_pipeline_from_config(pkg_salmon / "expected_import_project.yaml")
    bulk_manifest_path = Path(bulk_run_manifest["output_dir"]) / "run_manifest.json"
    assert bulk_manifest_path.is_file()

    # 2. Real scRNA 10x ingestion: run_scrna_validation reads 10x triad files
    _fastq(tmp_path / "star_real" / "fastq" / "S2_R1.fastq.gz")
    _fastq(tmp_path / "star_real" / "fastq" / "S2_R2.fastq.gz")
    sheet_star = _sheet(
        tmp_path / "star_real",
        header="sample\tfastq_1\tfastq_2",
        rows=["S1\tfastq/S1_R1.fastq.gz\tfastq/S1_R2.fastq.gz",
              "S2\tfastq/S2_R1.fastq.gz\tfastq/S2_R2.fastq.gz"],
    )
    pkg_star = _package(tmp_path / "star_real", "starsolo", samplesheet=sheet_star)
    sdirs_star, _ = _setup_valid_starsolo(pkg_star, samples=("S1", "S2"))

    # Execute real reading code via scrna_smoke.run_scrna_validation
    import ultimate.scrna_smoke as scrna_mod
    scrna_out = tmp_path / "scrna_downstream_run"
    star_run_manifest = scrna_mod.run_scrna_validation(
        input_path=sdirs_star[0],
        input_type="10x_mtx",
        output_dir=scrna_out,
        config={"modules": {"scrna": {"samples": [{"sample_id": "S1", "input_path": str(sdirs_star[0])},
                                                    {"sample_id": "S2", "input_path": str(sdirs_star[1])}]}}},
    )
    star_manifest_path = scrna_out / "run_manifest.json"
    assert star_manifest_path.is_file()

    # Both real manifests must be strictly accepted under dual-route verification
    res = runner.invoke(
        main,
        [
            "check-upstream-acceptance",
            "--package-dir", str(pkg_salmon), "--downstream-manifest", str(bulk_manifest_path),
            "--package-dir", str(pkg_star), "--downstream-manifest", str(star_manifest_path),
        ],
    )
    assert res.exit_code == 0, res.output
    payload = json.loads(res.output)
    assert payload["status"] == "pass"
    assert payload["failed_routes"] == []
    assert payload["missing_routes"] == []


def test_acceptance_counterexample_1_dimensions_mismatch(tmp_path: Path):
    """Counterexample 1: Upstream is 1x1, but downstream reports 71x88 -> dimension mismatch fails."""
    runner = CliRunner()
    pkg = _package(tmp_path, "salmon")
    matrix, m_valid = _setup_valid_salmon(pkg)
    # Upstream has 2 features x 1 sample. Tamper downstream to report 71 features x 88 samples.
    tampered_manifest = json.loads(m_valid.read_text(encoding="utf-8"))
    tampered_manifest["n_features"] = 71
    tampered_manifest["n_samples"] = 88
    tampered_manifest["raw_dimensions"] = {"n_features": 71, "n_samples": 88}
    m_bad = tmp_path / "m_bad_dims.json"
    m_bad.write_text(json.dumps(tampered_manifest), encoding="utf-8")

    res = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_bad), "--allow-single-route"],
    )
    assert res.exit_code != 0
    assert "downstream_dimensions_match_upstream" in res.output


def test_acceptance_counterexample_2_missing_feature_dimensions(tmp_path: Path):
    """Counterexample 2: Downstream only reports n_cells=1 without feature dimensions -> fails."""
    runner = CliRunner()
    pkg = _package(tmp_path, "starsolo")
    sdirs, m_valid = _setup_valid_starsolo(pkg, samples=("S1",))
    tampered_manifest = json.loads(m_valid.read_text(encoding="utf-8"))
    tampered_manifest["n_cells"] = 1
    tampered_manifest.pop("n_genes", None)
    tampered_manifest["raw_dimensions"] = {"n_cells": 1}  # feature dimension missing
    m_bad = tmp_path / "m_no_feature_dim.json"
    m_bad.write_text(json.dumps(tampered_manifest), encoding="utf-8")

    res = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_bad), "--allow-single-route"],
    )
    assert res.exit_code != 0
    assert "downstream_feature_dimension_present" in res.output


def test_acceptance_counterexample_3_hash_under_unrelated_filename(tmp_path: Path):
    """Counterexample 3: Correct SHA-256 bound to an unrelated filename must fail (no loose matching)."""
    runner = CliRunner()
    pkg = _package(tmp_path, "salmon")
    matrix, m_valid = _setup_valid_salmon(pkg)
    correct_sha = sha256_file(matrix)

    tampered_manifest = json.loads(m_valid.read_text(encoding="utf-8"))
    # Put correct hash under unrelated filename; counts_raw.tsv is unmapped
    tampered_manifest["input_checksums"] = {"completely_unrelated_dummy.txt": correct_sha}
    tampered_manifest.pop("input_checksum", None)
    m_bad = tmp_path / "m_unrelated_hash.json"
    m_bad.write_text(json.dumps(tampered_manifest), encoding="utf-8")

    res = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_bad), "--allow-single-route"],
    )
    assert res.exit_code != 0
    assert "checksum_not_bound_to_file:counts_raw.tsv" in res.output


def test_acceptance_counterexample_4_10x_triad_file_tampered(tmp_path: Path):
    """Counterexample 4: Any change to any of the 10x triad files (e.g. barcodes) fails checksum."""
    runner = CliRunner()
    pkg = _package(tmp_path, "starsolo")
    sdirs, m_valid = _setup_valid_starsolo(pkg, samples=("S1",))

    # Tamper with barcodes.tsv.gz after manifest was recorded
    barcodes_file = sdirs[0] / "barcodes.tsv.gz"
    with gzip.open(barcodes_file, "wt", encoding="utf-8") as h:
        h.write("TAMPERED_BARCODE_001\n")

    res = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_valid), "--allow-single-route"],
    )
    assert res.exit_code != 0
    assert "checksum_mismatch:barcodes.tsv.gz" in res.output


def test_m6_counterexample_missing_required_route(tmp_path: Path):
    """Single route fails final M6 acceptance under default require-both-routes."""
    pkg_salmon = _package(tmp_path / "s", "salmon")
    _, m_salmon = _setup_valid_salmon(pkg_salmon)
    runner = CliRunner()
    res = runner.invoke(
        main,
        ["check-upstream-acceptance", "--package-dir", str(pkg_salmon), "--downstream-manifest", str(m_salmon)],
    )
    assert res.exit_code != 0
    payload = json.loads(res.output)
    assert payload["status"] == "fail"
    assert "starsolo" in payload["missing_routes"]


def test_m6_counterexample_delivery_allowed_gate(tmp_path: Path):
    """delivery_allowed must be boolean False (not string, omitted, or True)."""
    runner = CliRunner()
    pkg = _package(tmp_path, "salmon")
    _, m_valid = _setup_valid_salmon(pkg)
    payload = json.loads(m_valid.read_text(encoding="utf-8"))

    for bad in (None, "false", True):
        p_bad = dict(payload)
        if bad is None:
            p_bad.pop("delivery_allowed", None)
        else:
            p_bad["delivery_allowed"] = bad
        m_file = tmp_path / f"m_gate_{bad}.json"
        m_file.write_text(json.dumps(p_bad), encoding="utf-8")
        res = runner.invoke(
            main,
            ["check-upstream-acceptance", "--package-dir", str(pkg), "--downstream-manifest", str(m_file), "--allow-single-route"],
        )
        assert res.exit_code != 0
        assert "delivery_allowed_boolean_false" in res.output


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
    assert 'grep -q $\'\\tfailed\\t\'' not in text
    assert 'not_ready=' in text
    assert "exit 1" in text
    assert "ULTIMATE_E2E_ALLOW_SKIP" in text
