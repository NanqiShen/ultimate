"""``ultimate doctor --module`` and the full handoff template scan.

``doctor`` must stay metadata-only: it may probe interpreters and import
specifications, but it must never run an analysis or write into a project.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from ultimate.cli import main
from ultimate.constants import MODULE_SPECS
from ultimate.diagnostics import inspect_environment, inspect_module
from ultimate.handoff_check import FORBIDDEN_TEMPLATE_TOKENS, run_handoff_check


# --- doctor --module ---------------------------------------------------------


def test_inspect_module_covers_known_modules():
    for module_name in ("rnaseq", "scrna", "proteomics", "clinical_assoc"):
        report = inspect_module(module_name)
        assert report["name"] == module_name
        assert report["status"] in {"pass", "fail"}
        assert report["checks"]


def test_inspect_module_rejects_unknown_module():
    report = inspect_module("not_a_module")
    assert report["status"] == "fail"
    assert report["reason"] == "unknown_module"
    assert "rnaseq" in report["known_modules"]


def test_inspect_module_checks_required_assets():
    report = inspect_module("rnaseq")
    names = {check["name"] for check in report["checks"]}
    assert "scripts:R" in names
    assert "slurm:ultimate_run.sbatch" in names
    assert "rscript" in names
    # The R backend must be declared for a module that uses DESeq2/edgeR.
    assert any(name.startswith("r:") for name in names)


def test_inspect_module_reports_missing_r_packages_as_warn_not_fail(tmp_path, monkeypatch):
    """A missing optional R package is a warning; it never fabricates success."""
    from ultimate import env_router

    fake = tmp_path / "Rscript"
    fake.write_text("#!/usr/bin/env sh\necho 'DESeq2\\tFALSE'\n", encoding="utf-8")
    fake.chmod(0o755)
    monkeypatch.setattr(env_router.shutil, "which", lambda name: str(fake) if name == "Rscript" else None)
    monkeypatch.delenv(env_router.ENV_RSCRIPT, raising=False)

    report = inspect_module("rnaseq")

    r_checks = [check for check in report["checks"] if check["name"].startswith("r:")]
    assert r_checks
    assert all(check["status"] in {"pass", "warn"} for check in r_checks)


def test_doctor_module_flag_is_wired(tmp_path):
    result = CliRunner().invoke(main, ["doctor", "--module", "rnaseq"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["module"]["name"] == "rnaseq"
    assert payload["scope"] == "metadata_only_not_production_approval"


def test_doctor_module_unknown_exits_nonzero():
    result = CliRunner().invoke(main, ["doctor", "--module", "definitely_not_a_module"])

    assert result.exit_code == 1
    payload = json.loads(result.output)
    assert payload["module"]["reason"] == "unknown_module"


def test_doctor_module_accepts_config(tmp_path):
    config = tmp_path / "project.yaml"
    config.write_text("project:\n  output_dir: ./runs\n  organism: human\n", encoding="utf-8")

    result = CliRunner().invoke(main, ["doctor", "--config", str(config), "--module", "scrna"])

    assert result.exit_code == 0, result.output


def test_doctor_module_does_not_create_output_dir(tmp_path, monkeypatch):
    config = tmp_path / "project.yaml"
    target = tmp_path / "uncreated_runs"
    config.write_text(f"project:\n  output_dir: {target}\n", encoding="utf-8")

    result = CliRunner().invoke(main, ["doctor", "--config", str(config), "--module", "rnaseq"])

    assert result.exit_code == 0, result.output
    assert not target.exists()


def test_inspect_environment_without_module_keeps_old_shape():
    report = inspect_environment()
    assert "module" not in report
    assert report["scope"] == "metadata_only_not_production_approval"


# --- handoff template scan ---------------------------------------------------


def test_forbidden_tokens_cover_pku_and_legacy_queues():
    for token in ("/share/home", "/shared/shen", "hpc-sbatch", "hpc-run", "queue = 'comp'", "queue = 'cpu'"):
        assert token in FORBIDDEN_TEMPLATE_TOKENS


def test_handoff_check_scans_every_template(tmp_path):
    manifest = run_handoff_check(root=tmp_path, output_dir=tmp_path / "hc")

    assert manifest["status"] == "ready", manifest["blockers"][:5]
    assert manifest["blocked_checks"] == 0
    assert len(manifest["scanned_handoff_templates"]) >= 30
    # The two nf-core handoffs keep their depth checks.
    assert {"nfcore_rnaseq", "nfcore_scrnaseq"}.issubset(set(manifest["checked_handoffs"]))


def test_handoff_check_reports_skipped_separately_from_blocked(tmp_path):
    manifest = run_handoff_check(root=tmp_path, output_dir=tmp_path / "hc")
    assert manifest["skipped_checks"] > 0
    assert all(":skipped" not in blocker for blocker in manifest["blockers"])


def test_handoff_check_blocks_a_reintroduced_pku_path(tmp_path, monkeypatch):
    """A template carrying /shared/shen again must be caught by the static scan."""
    from ultimate import handoff_check
    from ultimate import template_resources

    fake_root = tmp_path / "templates"
    handoff = fake_root / "handoffs" / "nfcore_rnaseq"
    handoff.mkdir(parents=True)
    (handoff / "README.md").write_text("# readme\n/shared/shen/2026/ultimate\n", encoding="utf-8")
    (handoff / "nextflow.config").write_text("process {\n  queue = 'comp'\n}\n", encoding="utf-8")
    monkeypatch.setattr(handoff_check, "template_root_candidates", lambda: [fake_root])

    manifest = run_handoff_check(root=tmp_path, output_dir=tmp_path / "hc")

    assert manifest["status"] == "blocked"
    joined = " ".join(manifest["blockers"])
    assert "forbidden literals" in joined
    assert "partition_from_env" in joined


def test_nfcore_nextflow_configs_are_partition_driven():
    from ultimate import env_router

    root = env_router.packaged_asset_candidates("templates", "handoffs")
    configs = [path / "nextflow.config" for path in root[0].glob("nfcore_*")]
    live = [path for path in configs if path.is_file()]
    assert live
    for config in live:
        text = config.read_text(encoding="utf-8")
        assert "ULTIMATE_SLURM_PARTITION" in text
        assert "queue = 'comp'" not in text
        assert "queue = 'cpu'" not in text


# --- thread budget on scrna commands ----------------------------------------


@pytest.mark.parametrize(
    "command",
    ["validate-scrna", "validate-scrna-pseudobulk-de", "annotate-scrna-celltypist"],
)
def test_scrna_commands_expose_threads(command: str):
    result = CliRunner().invoke(main, [command, "--help"])
    assert result.exit_code == 0
    assert "--threads" in result.output


def test_module_specs_used_by_doctor_exist():
    for module_name in ("rnaseq", "scrna"):
        assert module_name in MODULE_SPECS
