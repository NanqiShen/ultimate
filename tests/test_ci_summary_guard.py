"""The CI summary-line assertion must reject non-results.

This guard exists because a wrapper once reported success for a pytest usage
error. Any change that makes ``check_pytest_summary`` permissive reintroduces
that failure mode, so both acceptance and rejection paths are pinned.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "check_pytest_summary.py"

sys.path.insert(0, str(SCRIPT.parent))


def _load_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location("check_pytest_summary", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def summary_module():
    return _load_module()


def test_clean_run_passes(summary_module):
    text = "................................\n380 passed, 30 skipped, 149 warnings in 398.07s\n"
    assert summary_module.check(text) == []


def test_clean_run_parses_counts(summary_module):
    counts = summary_module.parse_summary("380 passed, 30 skipped in 398.07s")
    assert counts == {"passed": 380, "failed": 0, "error": 0, "skipped": 30, "found": 1}


def test_failed_tests_are_rejected(summary_module):
    problems = summary_module.check("1 failed, 379 passed in 12.00s")
    assert any("failed" in problem for problem in problems)


def test_errors_are_rejected(summary_module):
    problems = summary_module.check("2 errors, 380 passed in 12.00s")
    assert any("error" in problem for problem in problems)


def test_zero_passed_is_rejected(summary_module):
    problems = summary_module.check("0 passed, 5 skipped in 1.00s")
    assert any("no tests passed" in problem for problem in problems)


def test_usage_error_is_rejected_even_without_summary(summary_module):
    """The exact failure that was previously swallowed as exit_code=0."""
    text = (
        "ERROR: usage: python -m pytest [options] [file_or_dir] [file_or_dir] [...]\n"
        "python -m pytest: error: unrecognized arguments: --timeout=300\n"
    )
    problems = summary_module.check(text)
    assert any("unrecognized arguments" in problem for problem in problems)
    assert any("no pytest summary line" in problem for problem in problems)


def test_missing_summary_is_rejected(summary_module):
    problems = summary_module.check("just some noise without a result line\n")
    assert any("no pytest summary line" in problem for problem in problems)


def test_no_tests_ran_is_rejected(summary_module):
    problems = summary_module.check("no tests ran in 0.02s\n")
    assert any("no tests passed" in problem for problem in problems)


def test_cli_returns_zero_on_clean_log(tmp_path):
    log = tmp_path / "pytest.log"
    log.write_text("380 passed, 30 skipped in 398.07s\n", encoding="utf-8")
    completed = subprocess.run([sys.executable, str(SCRIPT), str(log)], capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "assertion: OK" in completed.stdout


def test_cli_returns_nonzero_on_usage_error(tmp_path):
    log = tmp_path / "pytest.log"
    log.write_text("error: unrecognized arguments: --timeout=300\n", encoding="utf-8")
    completed = subprocess.run([sys.executable, str(SCRIPT), str(log)], capture_output=True, text=True, check=False)
    assert completed.returncode == 1


def test_cli_returns_nonzero_on_missing_file(tmp_path):
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path / "absent.log")], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 2


def test_ci_workflow_uses_the_summary_script():
    workflow = (REPO_ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert "check_pytest_summary.py" in workflow
    # The workflow must not rely on the pytest exit code alone.
    assert "PYTHONPATH=src python -m pytest" in workflow
    assert "tests/test_environment_roots.py" not in workflow


def test_ci_matrix_covers_both_python_versions():
    workflow = (REPO_ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert '"3.11"' in workflow and '"3.12"' in workflow
    assert "doctor --module rnaseq" in workflow
    assert "doctor --module scrna" in workflow
    assert "handoff-check" in workflow
