"""Cross-environment routing: Rscript discovery and packaged asset lookup.

XieLab never exports conda/mamba/micromamba/Rscript on the default PATH, so the
four-step ``resolve_rscript`` priority is the only thing that lets an isolated
R 4.4 environment serve the limma/DESeq2/edgeR backends.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from ultimate import env_router
from ultimate.paths import ENV_ENV_ROOT, ENV_PLATFORM_ROOT, ENV_RUNTIME_ROOT


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for name in (
        ENV_PLATFORM_ROOT,
        ENV_RUNTIME_ROOT,
        ENV_ENV_ROOT,
        env_router.ENV_RSCRIPT,
        "ULTIMATE_ASSET_ROOT",
        "ULTIMATE_R_ENV_SCRNA",
        "ULTIMATE_R_ENV_RNASEQ",
    ):
        monkeypatch.delenv(name, raising=False)
    # PATH-based Rscript would mask the environment_root step.
    monkeypatch.setattr(env_router.shutil, "which", lambda name: None)


def _make_rscript(path: Path, executable: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/usr/bin/env sh\nexit 0\n", encoding="utf-8")
    path.chmod(0o755 if executable else 0o644)
    return path


# --- resolve_rscript: four-step priority -------------------------------------


def test_priority_1_module_config_wins(tmp_path, monkeypatch):
    monkeypatch.setenv(ENV_ENV_ROOT, str(tmp_path / "envs"))
    _make_rscript(tmp_path / "envs" / "ultimate-bio-r" / "bin" / "Rscript")
    configured = _make_rscript(tmp_path / "configured" / "Rscript")
    monkeypatch.setenv(env_router.ENV_RSCRIPT, str(_make_rscript(tmp_path / "env_override" / "Rscript")))

    assert env_router.resolve_rscript("scrna", {"rscript": str(configured)}) == configured.resolve()


def test_priority_1_module_config_accepts_de_backend_nesting(tmp_path):
    configured = _make_rscript(tmp_path / "de" / "Rscript")
    cfg = {"de_backend": {"rscript": str(configured)}}
    assert env_router.resolve_rscript("rnaseq", cfg) == configured.resolve()


def test_priority_1_non_executable_config_raises_而不是静默回退(tmp_path, monkeypatch):
    """An explicit but broken path must not be silently replaced by another R."""
    envs = tmp_path / "envs"
    monkeypatch.setenv(ENV_ENV_ROOT, str(envs))
    _make_rscript(envs / "ultimate-bio-r" / "bin" / "Rscript")
    not_executable = _make_rscript(tmp_path / "bad" / "Rscript", executable=False)

    with pytest.raises(ValueError, match="not an executable file"):
        env_router.resolve_rscript("scrna", {"rscript": str(not_executable)})

    missing = tmp_path / "absent" / "Rscript"
    with pytest.raises(ValueError, match="not an executable file"):
        env_router.resolve_rscript("scrna", {"rscript": str(missing)})

    described = env_router.describe_rscript("scrna", {"rscript": str(missing)})
    assert described["status"] == "missing"
    assert described["source"] == "module_config_invalid"


def test_priority_2_module_env_var_beats_global(tmp_path, monkeypatch):
    envs = tmp_path / "envs"
    monkeypatch.setenv(ENV_ENV_ROOT, str(envs))
    _make_rscript(envs / "ultimate-bio-r" / "bin" / "Rscript")
    module_r = _make_rscript(tmp_path / "module_r" / "Rscript")
    global_r = _make_rscript(tmp_path / "global_r" / "Rscript")
    monkeypatch.setenv("ULTIMATE_R_ENV_SCRNA", str(module_r))
    monkeypatch.setenv(env_router.ENV_RSCRIPT, str(global_r))

    assert env_router.resolve_rscript("scrna") == module_r.resolve()


def test_priority_2_global_env_var_used_when_module_var_absent(tmp_path, monkeypatch):
    monkeypatch.setenv(ENV_ENV_ROOT, str(tmp_path / "envs"))
    global_r = _make_rscript(tmp_path / "global_r" / "Rscript")
    monkeypatch.setenv(env_router.ENV_RSCRIPT, str(global_r))

    assert env_router.resolve_rscript("rnaseq") == global_r.resolve()


def test_priority_3_environment_root_preferred_envs_in_order(tmp_path, monkeypatch):
    envs = tmp_path / "envs"
    monkeypatch.setenv(ENV_ENV_ROOT, str(envs))
    rnaseq_r = _make_rscript(envs / "ultimate-rnaseq" / "bin" / "Rscript")
    # ultimate-bio-r is earlier in the default order but absent here.
    assert env_router.resolve_rscript("rnaseq") == rnaseq_r.resolve()


def test_priority_3_first_preferred_env_wins(tmp_path, monkeypatch):
    envs = tmp_path / "envs"
    monkeypatch.setenv(ENV_ENV_ROOT, str(envs))
    bio_r = _make_rscript(envs / "ultimate-bio-r" / "bin" / "Rscript")
    _make_rscript(envs / "ultimate-rnaseq" / "bin" / "Rscript")

    assert env_router.resolve_rscript("rnaseq") == bio_r.resolve()


def test_priority_4_path_fallback(tmp_path, monkeypatch):
    monkeypatch.setenv(ENV_ENV_ROOT, str(tmp_path / "envs"))
    on_path = _make_rscript(tmp_path / "path" / "Rscript")
    monkeypatch.setattr(env_router.shutil, "which", lambda name: str(on_path) if name == "Rscript" else None)

    assert env_router.resolve_rscript("scrna") == on_path.resolve()


def test_returns_none_when_nothing_found(tmp_path, monkeypatch):
    monkeypatch.setenv(ENV_ENV_ROOT, str(tmp_path / "envs"))
    assert env_router.resolve_rscript("scrna") is None
    described = env_router.describe_rscript("scrna")
    assert described["status"] == "missing"
    assert "not_found" == described["source"]


def test_module_env_variable_name_is_stable():
    assert env_router.module_env_variable("scrna") == "ULTIMATE_R_ENV_SCRNA"
    assert env_router.module_env_variable("cite_seq") == "ULTIMATE_R_ENV_CITE_SEQ"


def test_describe_rscript_reports_the_source(tmp_path, monkeypatch):
    envs = tmp_path / "envs"
    monkeypatch.setenv(ENV_ENV_ROOT, str(envs))
    _make_rscript(envs / "ultimate-bio-r" / "bin" / "Rscript")
    described = env_router.describe_rscript("scrna")
    assert described["status"] == "ready"
    assert described["source"] == "environment_root:ultimate-bio-r"

    described_cfg = env_router.describe_rscript("scrna", {"rscript": str(_make_rscript(tmp_path / "c" / "Rscript"))})
    assert described_cfg["source"] == "module_config"


# --- resolve_packaged_asset: source checkout and wheel ---------------------


def test_packaged_asset_found_in_source_checkout():
    long_path = env_router.resolve_packaged_asset("scripts", "R", "rnaseq_de_backend.R")
    assert long_path.is_file()
    assert long_path.name == "rnaseq_de_backend.R"


def test_packaged_asset_found_in_wheel_layout(tmp_path, monkeypatch):
    """A wheel install places data files under <prefix>/share/ultimate."""
    prefix = tmp_path / "venv"
    asset = prefix / "share" / "ultimate" / "scripts" / "R" / "wheel_only.R"
    asset.parent.mkdir(parents=True)
    asset.write_text("# bundled backend\n", encoding="utf-8")
    monkeypatch.setattr(env_router.sys, "prefix", str(prefix))

    resolved = env_router.resolve_packaged_asset("scripts", "R", "wheel_only.R")

    assert resolved == asset.resolve()


def test_ultimate_root_override_takes_precedence(tmp_path, monkeypatch):
    site = tmp_path / "site"
    asset = site / "scripts" / "R" / "site_only.R"
    asset.parent.mkdir(parents=True)
    asset.write_text("# site backend\n", encoding="utf-8")
    monkeypatch.setenv(ENV_PLATFORM_ROOT, str(site))

    assert env_router.resolve_packaged_asset("scripts", "R", "site_only.R") == asset.resolve()


def test_asset_root_override(tmp_path, monkeypatch):
    assets = tmp_path / "assets"
    asset = assets / "config" / "default.yaml"
    asset.parent.mkdir(parents=True)
    asset.write_text("project: {}\n", encoding="utf-8")
    monkeypatch.setenv("ULTIMATE_ASSET_ROOT", str(assets))

    assert env_router.resolve_packaged_asset("config", "default.yaml") == asset.resolve()


def test_missing_asset_raises_with_search_listing(tmp_path, monkeypatch):
    monkeypatch.setenv("ULTIMATE_ASSET_ROOT", str(tmp_path / "assets"))

    with pytest.raises(FileNotFoundError) as excinfo:
        env_router.resolve_packaged_asset("scripts", "R", "definitely_absent.R")

    assert "definitely_absent.R" in str(excinfo.value)
    assert "Searched:" in str(excinfo.value)
    assert env_router.find_packaged_asset("scripts", "R", "definitely_absent.R") is None


def test_packaged_asset_status_reports_candidates():
    status = env_router.packaged_asset_status("scripts", "R", "rnaseq_de_backend.R")
    assert status["status"] == "ready"
    assert status["path"].endswith("rnaseq_de_backend.R")
    assert status["searched"]


def test_convenience_accessors_point_at_real_files():
    assert env_router.scripts_r_dir() is not None and env_router.scripts_r_dir().is_dir()
    assert (env_router.scripts_r_dir() / "rnaseq_de_backend.R").is_file()
    assert env_router.default_config_path() is not None
    assert env_router.ultimate_run_sbatch_path() is not None


def test_every_bundled_r_backend_resolves():
    """Every scripts/R/*.R must be reachable through the router (wheel-safe)."""
    scripts_dir = env_router.scripts_r_dir()
    assert scripts_dir is not None
    for script in sorted(scripts_dir.glob("*.R")):
        assert env_router.resolve_packaged_asset("scripts", "R", script.name).is_file()
