"""CLI entry points required by the merged XieLab baseline.

The lineage merge had to keep two independently evolved CLI surfaces:

* the portable baseline (``doctor``, ``debug``, ``--threads``, ``--debug-log``),
  and
* the XieLab TFC line (``validate-scrna-pseudobulk-de``,
  ``annotate-scrna-celltypist``, native ``sbatch`` guidance).

These tests fail loudly if a future merge silently drops either side.
"""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from ultimate.cli import main


def test_portable_entries_are_registered() -> None:
    for command in ("doctor", "debug"):
        assert command in main.commands


def test_xielab_entries_are_registered() -> None:
    for command in ("validate-scrna-pseudobulk-de", "annotate-scrna-celltypist"):
        assert command in main.commands


def test_help_lists_both_lineages() -> None:
    result = CliRunner().invoke(main, ["--help"])
    assert result.exit_code == 0, result.output
    for command in ("doctor", "debug", "validate-scrna-pseudobulk-de", "annotate-scrna-celltypist"):
        assert command in result.output


def test_doctor_reports_metadata_only_scope(tmp_path: Path) -> None:
    config = tmp_path / "doctor.yaml"
    config.write_text("project:\n  output_dir: ./runs\n", encoding="utf-8")

    result = CliRunner().invoke(main, ["doctor", "--config", str(config)])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["scope"] == "metadata_only_not_production_approval"
    assert "checks" in payload and payload["checks"]


def test_doctor_fails_on_invalid_config(tmp_path: Path) -> None:
    config = tmp_path / "bad.yaml"
    config.write_text("project:\n  organism: not_an_organism\n", encoding="utf-8")

    result = CliRunner().invoke(main, ["doctor", "--config", str(config)])

    assert result.exit_code == 1
    assert json.loads(result.output)["status"] == "fail"


def test_debug_requires_config_or_trace() -> None:
    result = CliRunner().invoke(main, ["debug"])

    assert result.exit_code != 0
    assert "provide --config and/or --trace" in result.output


def test_debug_summarizes_trace_without_replaying(tmp_path: Path) -> None:
    trace = tmp_path / "trace.jsonl"
    trace.write_text(
        "\n".join(
            json.dumps(event)
            for event in (
                {"stage": "config", "status": "started"},
                {"stage": "config", "status": "completed"},
                {"stage": "run", "status": "failed"},
            )
        )
        + "\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(main, ["debug", "--trace", str(trace)])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["trace"]["first_failed_stage"] == "run"
    assert payload["trace"]["trace_events"]["failed"] == 1


def test_debug_rejects_unknown_trace_event(tmp_path: Path) -> None:
    trace = tmp_path / "trace.jsonl"
    trace.write_text('{"stage": "unknown", "status": "started"}\n', encoding="utf-8")

    result = CliRunner().invoke(main, ["debug", "--trace", str(trace)])

    assert result.exit_code != 0
    assert "unknown trace event" in result.output


def test_validate_scrna_pseudobulk_de_requires_input(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        main,
        [
            "validate-scrna-pseudobulk-de",
            "--output-dir",
            str(tmp_path / "out"),
            "--cell-type",
            "T",
            "--control",
            "ctrl",
            "--case",
            "stim",
        ],
    )

    assert result.exit_code != 0
    # The command must explain that an input source is required rather than
    # producing an empty analysis.
    assert "sce" in result.output.lower() or "counts" in result.output.lower()


def test_validate_scrna_pseudobulk_de_reports_missing_input_file(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        main,
        [
            "validate-scrna-pseudobulk-de",
            "--output-dir",
            str(tmp_path / "out"),
            "--counts",
            str(tmp_path / "absent.tsv"),
            "--cell-type",
            "T",
            "--control",
            "ctrl",
            "--case",
            "stim",
        ],
    )

    assert result.exit_code != 0


def test_annotate_scrna_celltypist_requires_model_and_input(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        main,
        [
            "annotate-scrna-celltypist",
            "--input-h5ad",
            str(tmp_path / "missing.h5ad"),
            "--output-dir",
            str(tmp_path / "out"),
            "--model",
            str(tmp_path / "missing.pkl"),
        ],
    )

    assert result.exit_code != 0


def test_run_accepts_threads_and_debug_log_options() -> None:
    help_result = CliRunner().invoke(main, ["run", "--help"])

    assert help_result.exit_code == 0
    assert "--threads" in help_result.output
    assert "--debug-log" in help_result.output


def test_debug_log_refuses_to_overwrite(tmp_path: Path) -> None:
    """``--debug-log`` opens in exclusive mode so a previous trace is preserved."""
    config = tmp_path / "project.yaml"
    config.write_text("project:\n  organism: not_an_organism\n", encoding="utf-8")
    log = tmp_path / "trace.jsonl"
    log.write_text("existing\n", encoding="utf-8")

    result = CliRunner().invoke(
        main, ["run", "--config", str(config), "--debug-log", str(log)]
    )

    assert result.exit_code != 0
    assert log.read_text(encoding="utf-8") == "existing\n"
