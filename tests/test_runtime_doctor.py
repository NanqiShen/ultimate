import io
import json
import os
from pathlib import Path

import pytest
from click.testing import CliRunner

from ultimate.cli import main
from ultimate.runtime import RunTrace, apply_thread_budget, thread_budget


def test_budget_never_exceeds_slurm_or_affinity(monkeypatch):
    monkeypatch.setattr(os, "cpu_count", lambda: 96)
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: {0, 1, 2}, raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("ULTIMATE_THREADS", "7")
    cfg = {"runtime": {"threads": 5}}
    result = thread_budget(cfg, 8)
    assert (result.threads, result.requested, result.source, result.capped) == (2, 8, "--threads", True)
    assert thread_budget(cfg).source == "ULTIMATE_THREADS"
    monkeypatch.delenv("ULTIMATE_THREADS")
    assert thread_budget(cfg).source == "runtime.threads"
    monkeypatch.delenv("SLURM_CPUS_PER_TASK")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    assert thread_budget({}).threads == 1


@pytest.mark.parametrize("value", [-2, "2.0", True])
def test_invalid_thread_setting_is_rejected(value):
    with pytest.raises(ValueError, match="positive integer"):
        thread_budget({}, value)


def test_thread_env_is_restored_after_failure(monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "12")
    monkeypatch.delenv("RCPP_PARALLEL_NUM_THREADS", raising=False)
    with pytest.raises(RuntimeError):
        with apply_thread_budget(thread_budget({}, 1)):
            assert os.environ["OMP_NUM_THREADS"] == "1"
            assert os.environ["RCPP_PARALLEL_NUM_THREADS"] == "1"
            raise RuntimeError("synthetic failure")
    assert os.environ["OMP_NUM_THREADS"] == "12"
    assert "RCPP_PARALLEL_NUM_THREADS" not in os.environ


def test_doctor_and_debug_are_read_only(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    path = tmp_path / "project.yaml"
    path.write_text("project:\n  output_dir: ./uncreated\n")
    r = CliRunner().invoke(main, ["doctor", "--config", str(path), "--threads", "8"])
    assert r.exit_code == 0, r.output
    assert json.loads(r.output)["threads"]["threads"] <= 2
    assert not (tmp_path / "uncreated").exists()
    path.write_text("project: [\n")
    r = CliRunner().invoke(main, ["doctor", "--config", str(path)])
    assert r.exit_code == 1 and json.loads(r.output)["status"] == "fail"


def test_trace_sections_and_existing_log(tmp_path):
    config = tmp_path / "invalid.yaml"
    config.write_text("project: [\n")
    out = tmp_path / "trace.jsonl"
    r = CliRunner().invoke(main, ["run", "--config", str(config), "--debug-log", str(out)])
    assert r.exit_code != 0
    r = CliRunner().invoke(main, ["debug", "--trace", str(out)])
    assert json.loads(r.output)["trace"]["first_failed_stage"] == "config"
    original = out.read_bytes()
    r = CliRunner().invoke(main, ["run", "--config", str(config), "--debug-log", str(out)])
    assert r.exit_code != 0 and out.read_bytes() == original
    stream = io.StringIO()
    with pytest.raises(RuntimeError):
        with RunTrace(stream).stage("module", module="test"):
            raise RuntimeError("/private/path should not appear")
    assert "/private/path" not in stream.getvalue()
