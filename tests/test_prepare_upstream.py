"""M5 upstream packaging: `ultimate prepare-upstream` job-package generation.

Pins the guarantees the operator relies on before submitting anything to Slurm:
the package is generated (never auto-submitted), FASTQ inputs stay read-only,
threads are bound to the run budget, and a missing index / FASTQ / samplesheet
blocks with a concrete reason instead of producing a broken script.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from ultimate.cli import main
from ultimate.prepare_upstream import (
    INDEX_MARKERS,
    ROUTES,
    _read_samplesheet,
    prepare_upstream_package,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------


def _make_fastq(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x1f\x8b\x08\x00fake-fastq")
    return path


def _make_salmon_index(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "info.json").write_text('{"salmon_version": "1.10.0"}', encoding="utf-8")
    return root


def _make_star_index(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "SAindex").write_text("binary", encoding="utf-8")
    (root / "Genome").write_text("binary", encoding="utf-8")
    return root


def _samplesheet(tmp_path: Path, *, paired: bool = True, samples=("S1", "S2")) -> Path:
    fastq_dir = tmp_path / "fastq"
    rows = ["sample\tfastq_1\tfastq_2" if paired else "sample\tfastq_1"]
    for sample in samples:
        _make_fastq(fastq_dir / f"{sample}_R1.fastq.gz")
        if paired:
            _make_fastq(fastq_dir / f"{sample}_R2.fastq.gz")
            rows.append(f"{sample}\tfastq/{sample}_R1.fastq.gz\tfastq/{sample}_R2.fastq.gz")
        else:
            rows.append(f"{sample}\tfastq/{sample}_R1.fastq.gz")
    path = tmp_path / "samplesheet.tsv"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# route coverage: salmon + starsolo packages
# --------------------------------------------------------------------------


def test_route_table_covers_salmon_and_starsolo():
    assert set(ROUTES) == {"salmon", "starsolo"}
    assert ROUTES["salmon"]["module"] == "rnaseq"
    assert ROUTES["starsolo"]["module"] == "scrna"
    assert set(INDEX_MARKERS) == {"salmon", "starsolo"}


@pytest.mark.parametrize(("route", "index_factory"), [("salmon", _make_salmon_index), ("starsolo", _make_star_index)])
def test_prepare_upstream_generates_complete_package(tmp_path: Path, route: str, index_factory):
    samplesheet = _samplesheet(tmp_path)
    index = index_factory(tmp_path / "index")
    out = tmp_path / "pkg"

    result = prepare_upstream_package(
        route=route, samplesheet=samplesheet, output_dir=out, index=index, threads=2,
    )

    assert result.status == "ready", result.blocked_reason
    assert result.threads == 2
    assert result.n_samples == 2
    for key in ("run_upstream_sbatch", "expected_import_project", "checksum_manifest", "job_package_manifest"):
        assert Path(result.artifacts[key]).is_file(), key

    script = Path(result.artifacts["run_upstream_sbatch"]).read_text(encoding="utf-8")
    # Slurm directives must sit in the leading comment block, before code.
    body = script.splitlines()
    first_code = next(i for i, line in enumerate(body) if line.strip() and not line.startswith("#"))
    sbatch_indexes = [i for i, line in enumerate(body) if line.startswith("#SBATCH")]
    assert sbatch_indexes and max(sbatch_indexes) < first_code
    assert "#SBATCH --account=xielab" in script
    assert "#SBATCH --partition=compute_pro" in script
    assert "#SBATCH --qos=normal" in script
    # Thread binding is consistent with the requested budget.
    assert "#SBATCH --cpus-per-task=2" in script
    assert "THREADS=${SLURM_CPUS_PER_TASK:-2}" in script
    assert "THREADS=$((SLURM_CPUS_PER_TASK" not in script

    manifest = json.loads(Path(result.artifacts["job_package_manifest"]).read_text(encoding="utf-8"))
    assert manifest["auto_submitted"] is False
    assert manifest["raw_read_only"] is True
    assert manifest["route"] == route
    assert manifest["submit_command"].startswith("sbatch ")


def test_generated_script_passes_bash_syntax(tmp_path: Path):
    import subprocess

    samplesheet = _samplesheet(tmp_path)
    index = _make_salmon_index(tmp_path / "index")
    result = prepare_upstream_package(
        route="salmon", samplesheet=samplesheet, output_dir=tmp_path / "pkg", index=index, threads=1,
    )
    script = Path(result.artifacts["run_upstream_sbatch"])
    completed = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    assert script.stat().st_mode & 0o111  # executable bit set for review


def test_expected_import_project_targets_route_output(tmp_path: Path):
    samplesheet = _samplesheet(tmp_path, paired=False)
    result = prepare_upstream_package(
        route="salmon", samplesheet=samplesheet, output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"),
    )
    text = Path(result.artifacts["expected_import_project"]).read_text(encoding="utf-8")
    assert "rnaseq:" in text and "counts_raw.tsv" in text

    result2 = prepare_upstream_package(
        route="starsolo", samplesheet=samplesheet, output_dir=tmp_path / "pkg2",
        index=_make_star_index(tmp_path / "index2"),
    )
    text2 = Path(result2.artifacts["expected_import_project"]).read_text(encoding="utf-8")
    assert "scrna:" in text2 and "10x_mtx" in text2


# --------------------------------------------------------------------------
# read-only input validation
# --------------------------------------------------------------------------


def test_checksum_manifest_records_inputs_as_read_only(tmp_path: Path):
    samplesheet = _samplesheet(tmp_path)
    import csv

    fastq = tmp_path / "fastq" / "S1_R1.fastq.gz"
    original = fastq.read_bytes()
    result = prepare_upstream_package(
        route="salmon", samplesheet=samplesheet, output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"),
    )
    rows = list(csv.DictReader(Path(result.artifacts["checksum_manifest"]).open(encoding="utf-8"), delimiter="\t"))
    assert rows and all(row["read_only"] == "true" for row in rows)
    assert {row["role"] for row in rows} == {"fastq_1", "fastq_2"}
    assert all(len(row["sha256"]) == 64 for row in rows)
    # The input file must be byte-identical after packaging.
    assert fastq.read_bytes() == original


def test_missing_fastq_blocks(tmp_path: Path):
    samplesheet = _samplesheet(tmp_path)
    (tmp_path / "fastq" / "S2_R1.fastq.gz").unlink()
    result = prepare_upstream_package(
        route="salmon", samplesheet=samplesheet, output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"),
    )
    assert result.status == "blocked"
    assert "missing_fastq:S2:fastq_1" in result.blocked_reason
    assert not (tmp_path / "pkg" / "run_upstream.sbatch").exists()


def test_non_fastq_input_blocks(tmp_path: Path):
    samplesheet = tmp_path / "samplesheet.tsv"
    (tmp_path / "fastq").mkdir()
    (tmp_path / "fastq" / "S1_R1.bam").write_text("not-fastq", encoding="utf-8")
    samplesheet.write_text("sample\tfastq_1\nS1\tfastq/S1_R1.bam\n", encoding="utf-8")
    result = prepare_upstream_package(
        route="salmon", samplesheet=samplesheet, output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"),
    )
    assert result.status == "blocked"
    assert "not_a_fastq_file" in result.blocked_reason


def test_missing_samplesheet_blocks(tmp_path: Path):
    result = prepare_upstream_package(
        route="salmon", samplesheet=tmp_path / "absent.tsv", output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"),
    )
    assert result.status == "blocked"
    assert "missing_samplesheet" in result.blocked_reason


def test_empty_samplesheet_blocks(tmp_path: Path):
    path = tmp_path / "s.tsv"
    path.write_text("sample\tfastq_1\n", encoding="utf-8")
    result = prepare_upstream_package(
        route="salmon", samplesheet=path, output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"),
    )
    assert result.status == "blocked"
    assert "empty_samplesheet" in result.blocked_reason


# --------------------------------------------------------------------------
# index requirements
# --------------------------------------------------------------------------


def test_missing_index_blocks(tmp_path: Path):
    result = prepare_upstream_package(
        route="salmon", samplesheet=_samplesheet(tmp_path), output_dir=tmp_path / "pkg", index=None,
    )
    assert result.status == "blocked"
    assert "missing_index" in result.blocked_reason


def test_index_directory_without_markers_blocks(tmp_path: Path):
    empty_index = tmp_path / "not_an_index"
    empty_index.mkdir()
    result = prepare_upstream_package(
        route="salmon", samplesheet=_samplesheet(tmp_path), output_dir=tmp_path / "pkg", index=empty_index,
    )
    assert result.status == "blocked"
    assert "index_not_usable_for_salmon" in result.blocked_reason


def test_star_route_rejects_salmon_index(tmp_path: Path):
    result = prepare_upstream_package(
        route="starsolo", samplesheet=_samplesheet(tmp_path, paired=False),
        output_dir=tmp_path / "pkg", index=_make_salmon_index(tmp_path / "sidx"),
    )
    assert result.status == "blocked"
    assert "index_not_usable_for_starsolo" in result.blocked_reason


def test_unsupported_route_blocks(tmp_path: Path):
    result = prepare_upstream_package(
        route="nextflow", samplesheet=_samplesheet(tmp_path), output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"),
    )
    assert result.status == "blocked"
    assert "unsupported_route" in result.blocked_reason


# --------------------------------------------------------------------------
# samplesheet parsing
# --------------------------------------------------------------------------


def test_samplesheet_requires_sample_and_fastq_columns(tmp_path: Path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("id\tfastq_1\nS1\tfastq/S1_R1.fastq.gz\n", encoding="utf-8")
    rows, error = _read_samplesheet(bad)
    assert rows == []
    assert "samplesheet_missing_sample_at_row" in error


def test_samplesheet_accepts_sample_id_alias(tmp_path: Path):
    path = tmp_path / "s.tsv"
    path.write_text("sample_id\tfq1\nS1\tfastq/S1_R1.fastq.gz\n", encoding="utf-8")
    rows, error = _read_samplesheet(path)
    assert error == ""
    assert rows[0]["sample_id"] == "S1"


# --------------------------------------------------------------------------
# thread budget
# --------------------------------------------------------------------------


def test_threads_are_capped_by_budget(tmp_path: Path, monkeypatch):
    import os

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.delenv("ULTIMATE_THREADS", raising=False)
    samplesheet = _samplesheet(tmp_path)
    result = prepare_upstream_package(
        route="salmon", samplesheet=samplesheet, output_dir=tmp_path / "pkg",
        index=_make_salmon_index(tmp_path / "index"), threads=8,
    )
    assert result.status == "ready"
    assert result.threads <= 2


# --------------------------------------------------------------------------
# CLI wiring
# --------------------------------------------------------------------------


def test_cli_prepare_upstream_ready(tmp_path: Path):
    samplesheet = _samplesheet(tmp_path)
    index = _make_salmon_index(tmp_path / "index")
    result = CliRunner().invoke(
        main,
        ["prepare-upstream", "--route", "salmon", "--samplesheet", str(samplesheet),
         "--index", str(index), "--output-dir", str(tmp_path / "pkg")],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["status"] == "ready"
    assert payload["route"] == "salmon"


def test_cli_prepare_upstream_blocked_exits_nonzero(tmp_path: Path):
    result = CliRunner().invoke(
        main,
        ["prepare-upstream", "--route", "salmon", "--samplesheet", str(_samplesheet(tmp_path)),
         "--output-dir", str(tmp_path / "pkg")],
    )
    assert result.exit_code != 0
    assert "missing_index" in result.output


# --------------------------------------------------------------------------
# triage wiring + repo hygiene
# --------------------------------------------------------------------------


def test_triage_recommends_prepare_upstream_for_fastq():
    from ultimate.triage import _prepare_upstream_hints

    hints = _prepare_upstream_hints(["rnaseq"])
    assert hints and hints[0]["route"] == "salmon"
    assert "prepare-upstream" in hints[0]["message"]
    assert hints[0]["status"] in {"available", "needs_dependency"}

    hints_scrna = _prepare_upstream_hints(["scrna"])
    assert hints_scrna[0]["route"] == "starsolo"


def test_new_slurm_scripts_have_xielab_directives():
    scripts = [
        "slurm/xielab_setup_upstream_env.sbatch",
        "slurm/upstream_bulk_rnaseq_salmon.sbatch",
        "slurm/upstream_scrna_starsolo.sbatch",
        "slurm/xielab_e2e_regression_suite.sbatch",
    ]
    for name in scripts:
        text = (REPO_ROOT / name).read_text(encoding="utf-8")
        assert "#SBATCH --account=xielab" in text, name
        assert "#SBATCH --partition=compute_pro" in text, name
        assert "#SBATCH --qos=normal" in text, name
        assert "hpc-sbatch" not in text and "hpc-run" not in text, name
        assert "SLURM_CPUS_PER_TASK" in text, name


def test_e2e_suite_never_auto_submits():
    text = (REPO_ROOT / "slurm" / "xielab_e2e_regression_suite.sbatch").read_text(encoding="utf-8")
    code_lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    assert not [line for line in code_lines if line.strip().startswith("sbatch ")]


def test_upstream_env_definition_covers_required_tools():
    import yaml

    payload = yaml.safe_load((REPO_ROOT / "envs" / "xielab.upstream.yml").read_text(encoding="utf-8"))
    deps = [str(item) for item in payload["dependencies"]]
    for tool in ("fastqc", "multiqc", "fastp", "salmon", "star", "subread", "samtools", "alevin-fry"):
        assert tool in deps, tool


def test_operations_doc_exists_and_documents_key_facts():
    text = (REPO_ROOT / "docs" / "XIELAB_OPERATIONS.md").read_text(encoding="utf-8")
    for fact in ("ultimate_runtime", "ultimate-upstream", "compute_pro", "xielab", "500 GB", "/dshare/rawdata"):
        assert fact in text, fact
    # PKU paths must never appear; mentioning that hpc-sbatch does not exist is
    # the intended negative statement and is allowed.
    assert "/shared/shen" not in text
    assert "hpc-sbatch" not in text.replace("集群**没有** `hpc-sbatch` / `hpc-run`", "")
