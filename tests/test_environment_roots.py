"""Environment roots are resolved through ``ultimate.paths`` everywhere.

These tests pin the three call sites that used to hard-code ``.conda/envs``:
the preflight dependency probes, the tool audit and the single-cell capability
audit. Without them a future refactor could silently reintroduce a
checkout-relative layout that does not exist on the cluster.
"""

from __future__ import annotations

from pathlib import Path

from ultimate import paths
from ultimate.preflight import _env_bin_roots, _env_root
from ultimate.singlecell_audit import CAPABILITIES, _env_paths
from ultimate.tool_registry import _env_paths as _tool_env_paths


def _clear(monkeypatch):
    for name in (paths.ENV_PLATFORM_ROOT, paths.ENV_RUNTIME_ROOT, paths.ENV_ENV_ROOT, paths.ENV_CACHE_ROOT):
        monkeypatch.delenv(name, raising=False)


def test_preflight_env_root_sits_in_runtime_sibling(tmp_path, monkeypatch):
    _clear(monkeypatch)
    config = {"project": {"server_root": str(tmp_path / "ultimate")}}
    assert _env_root(config) == (tmp_path / "ultimate_runtime" / "envs").resolve()
    bins = _env_bin_roots(config)
    assert bins and all(path.is_absolute() for path in bins)
    assert bins[0] == (tmp_path / "ultimate_runtime" / "envs" / "ultimate-core" / "bin")


def test_preflight_env_root_independent_of_project_noise(tmp_path, monkeypatch):
    """A portable checkout whose directory is not named ``ultimate`` nests the runtime."""
    _clear(monkeypatch)
    config = {"project": {"server_root": str(tmp_path / "portable-checkout")}}
    assert _env_root(config) == (tmp_path / "portable-checkout" / "ultimate_runtime" / "envs").resolve()


def test_preflight_env_root_honours_override(tmp_path, monkeypatch):
    _clear(monkeypatch)
    override = tmp_path / "managed envs"
    monkeypatch.setenv(paths.ENV_ENV_ROOT, str(override))
    config = {"project": {"server_root": str(tmp_path / "ultimate")}}
    assert _env_root(config) == override.resolve()


def test_preflight_env_root_without_server_root(monkeypatch):
    _clear(monkeypatch)
    assert _env_root({}) == paths.environment_root()


def test_tool_registry_env_paths_use_environment_root(tmp_path, monkeypatch):
    _clear(monkeypatch)
    root = tmp_path / "ultimate"
    root.mkdir()
    resolved = _tool_env_paths(root)
    assert resolved
    expected_root = (tmp_path / "ultimate_runtime" / "envs").resolve()
    for env, path in resolved.items():
        assert path == expected_root / env


def test_tool_registry_respects_env_override(tmp_path, monkeypatch):
    _clear(monkeypatch)
    override = tmp_path / "shared envs"
    monkeypatch.setenv(paths.ENV_ENV_ROOT, str(override))
    resolved = _tool_env_paths(tmp_path / "ultimate")
    assert all(path == override.resolve() / env for env, path in resolved.items())


def test_singlecell_audit_env_paths_translate_placeholders(tmp_path, monkeypatch):
    _clear(monkeypatch)
    root = tmp_path / "ultimate"
    root.mkdir()
    env_root = (tmp_path / "ultimate_runtime" / "envs").resolve()
    capability = next(cap for cap in CAPABILITIES if cap.key == "scatac")
    resolved = _env_paths(root, capability)
    assert resolved[0] == env_root / capability.env_name
    # Alternate environment paths must resolve under the same root, not the code root.
    assert (env_root / "ultimate-scatac-py") in resolved
    assert (env_root / "ultimate-scatac-r") in resolved
    assert not any("/.conda/" in str(path) for path in resolved)


def test_singlecell_audit_external_workbench_paths_stay_out_of_tree(tmp_path, monkeypatch):
    """The external singlecell_workbench keeps its own conda envs by design."""
    _clear(monkeypatch)
    external = tmp_path / "external"
    monkeypatch.setenv(paths.ENV_EXTERNAL_ROOT, str(external))
    capability = next(cap for cap in CAPABILITIES if cap.key == "scrna")
    resolved = _env_paths(tmp_path / "ultimate", capability)
    assert (external.resolve() / "singlecell_workbench" / ".conda" / "envs" / "scw-py311") in resolved


def test_every_capability_resolves_to_absolute_paths(tmp_path, monkeypatch):
    _clear(monkeypatch)
    root = tmp_path / "ultimate"
    root.mkdir()
    for capability in CAPABILITIES:
        for path in _env_paths(root, capability):
            assert isinstance(path, Path)
            assert path.is_absolute(), f"{capability.key} produced a relative env path: {path}"
            assert "{" not in str(path), f"{capability.key} left an unresolved placeholder: {path}"
