"""Unified path resolution: checkout layout plus XieLab runtime layout."""

from pathlib import Path

import pytest

from ultimate import paths


@pytest.fixture(autouse=True)
def _clean_path_env(monkeypatch):
    for name in (
        paths.ENV_PLATFORM_ROOT,
        paths.ENV_RESOURCES,
        paths.ENV_EXTERNAL_ROOT,
        paths.ENV_RUNTIME_ROOT,
        paths.ENV_CACHE_ROOT,
        paths.ENV_ENV_ROOT,
        paths.ENV_REFERENCE_ROOT,
    ):
        monkeypatch.delenv(name, raising=False)


def test_project_root_follows_source_checkout():
    root = paths.project_root()
    assert (root / "config/default.yaml").is_file()
    assert (root / "pyproject.toml").is_file()


def test_project_root_prefers_explicit_override(tmp_path, monkeypatch):
    override = tmp_path / "site root"
    monkeypatch.setenv(paths.ENV_PLATFORM_ROOT, str(override))
    assert paths.project_root() == override.resolve()
    # The compatibility alias must never drift from the primary resolver.
    assert paths.default_platform_root() == paths.project_root()


def test_resource_and_external_roots_default_under_project_root():
    root = paths.project_root()
    assert paths.resource_root() == (root / "resources").resolve()
    assert paths.external_root() == (root / "external").resolve()


@pytest.mark.parametrize(
    ("variable", "resolver", "suffix"),
    [
        (paths.ENV_RESOURCES, "resource_root", "resources"),
        (paths.ENV_EXTERNAL_ROOT, "external_root", "external"),
        (paths.ENV_RUNTIME_ROOT, "runtime_root", "runtime"),
        (paths.ENV_ENV_ROOT, "environment_root", "envs"),
        (paths.ENV_CACHE_ROOT, "cache_root", "cache"),
    ],
)
def test_each_root_honours_its_environment_variable(tmp_path, monkeypatch, variable, resolver, suffix):
    override = tmp_path / "override directory" / suffix
    monkeypatch.setenv(variable, str(override))
    resolved = getattr(paths, resolver)()
    assert resolved == override.resolve()
    assert resolved.is_absolute()


def test_runtime_root_sits_beside_classic_ultimate_checkout(tmp_path):
    project = tmp_path / "site" / "ultimate"
    project.mkdir(parents=True)
    assert paths.runtime_root(project) == (tmp_path / "site" / "ultimate_runtime").resolve()


@pytest.mark.parametrize("name", ["ultimate-code", "analysis", "scrna"])
def test_runtime_root_nests_inside_non_canonical_checkout(tmp_path, name):
    project = tmp_path / name
    project.mkdir(parents=True)
    assert paths.runtime_root(project) == (project / "ultimate_runtime").resolve()


def test_environment_and_cache_roots_derive_from_runtime_root(tmp_path, monkeypatch):
    runtime = tmp_path / "runtime"
    monkeypatch.setenv(paths.ENV_RUNTIME_ROOT, str(runtime))
    assert paths.environment_root() == (runtime / "envs").resolve()
    assert paths.cache_root() == (runtime / "cache").resolve()


def test_environment_root_independent_override_wins_over_runtime_root(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ENV_RUNTIME_ROOT, str(tmp_path / "runtime"))
    envs = tmp_path / "shared envs"
    monkeypatch.setenv(paths.ENV_ENV_ROOT, str(envs))
    assert paths.environment_root() == envs.resolve()
    assert paths.default_environment_root() == envs.resolve()
    # cache_root still follows the runtime root.
    assert paths.cache_root() == (tmp_path / "runtime" / "cache").resolve()


def test_xielab_layout_resolves_from_scratch(tmp_path, monkeypatch):
    """End-to-end check of the documented XieLab deployment shape."""
    code = tmp_path / "site" / "ultimate"
    code.mkdir(parents=True)
    monkeypatch.setenv(paths.ENV_PLATFORM_ROOT, str(code))
    assert paths.project_root() == code.resolve()
    assert paths.runtime_root() == (tmp_path / "site" / "ultimate_runtime").resolve()
    assert paths.environment_root().name == "envs"
    assert paths.cache_root().name == "cache"


def test_env_value_with_shell_variable_is_expanded(tmp_path, monkeypatch):
    monkeypatch.setenv("ULTIMATE_TEST_BASE", str(tmp_path))
    monkeypatch.setenv(paths.ENV_CACHE_ROOT, "$ULTIMATE_TEST_BASE/cache")
    assert paths.cache_root() == (tmp_path / "cache").resolve()


def test_blank_environment_value_is_ignored(tmp_path, monkeypatch):
    """An empty variable is a common shell accident and must not create ``Path('')``."""
    monkeypatch.setenv(paths.ENV_RUNTIME_ROOT, "")
    assert paths.runtime_root(tmp_path) == (tmp_path / "ultimate_runtime").resolve()


def test_resolve_path_handles_absolute_relative_and_variables(tmp_path, monkeypatch):
    monkeypatch.setenv("ULTIMATE_TEST_BASE", str(tmp_path))
    assert paths.resolve_path(tmp_path, "sub") == (tmp_path / "sub").resolve()
    assert paths.resolve_path(tmp_path, tmp_path / "abs") == (tmp_path / "abs").resolve()
    assert paths.resolve_path(tmp_path, "$ULTIMATE_TEST_BASE/var") == (tmp_path / "var").resolve()
    assert paths.expand_path("~/data").is_absolute()
