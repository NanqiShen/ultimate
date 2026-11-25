"""Deployment paths, resolved at runtime rather than pinned to a cluster.

The module is the single place where Ultimate maps environment overrides onto
filesystem locations. It merges two previously divergent conventions:

* the portable checkout layout (``project_root`` / ``resource_root`` /
  ``external_root``) used by source checkouts and site overrides, and
* the XieLab runtime layout (``runtime_root`` / ``environment_root`` /
  ``cache_root``) that keeps environments, caches, logs and manifests outside
  the read-only code root.

Both families are kept side by side so callers never need to guess. Every
resolver honours an environment variable first and only then falls back to a
deterministic default derived from :func:`project_root`.
"""

from __future__ import annotations

import os
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[2]

ENV_PLATFORM_ROOT = "ULTIMATE_ROOT"
ENV_RESOURCES = "ULTIMATE_RESOURCES"
ENV_EXTERNAL_ROOT = "ULTIMATE_EXTERNAL_ROOT"
ENV_RUNTIME_ROOT = "ULTIMATE_RUNTIME_ROOT"
ENV_CACHE_ROOT = "ULTIMATE_CACHE_ROOT"
ENV_ENV_ROOT = "ULTIMATE_ENV_ROOT"
ENV_REFERENCE_ROOT = "ULTIMATE_REFERENCE_ROOT"

RUNTIME_DIRNAME = "ultimate_runtime"


def expand_path(value: str | Path) -> Path:
    """Expand ``~`` and environment variables without touching the filesystem."""
    return Path(os.path.expandvars(str(value))).expanduser()


def resolve_path(base_dir: Path, value: str | Path) -> Path:
    """Resolve ``value`` against ``base_dir`` when it is relative."""
    candidate = expand_path(value)
    return candidate.resolve() if candidate.is_absolute() else (base_dir / candidate).resolve()


def _env_path(name: str) -> Path | None:
    value = os.environ.get(name)
    if not value:
        return None
    return expand_path(value).resolve()


def project_root() -> Path:
    """Execution root. Supports source checkouts and explicit site overrides."""
    value = _env_path(ENV_PLATFORM_ROOT)
    if value is not None:
        return value
    # Installed wheels live under site-packages, not a project checkout.
    if (CODE_ROOT / "pyproject.toml").is_file() and (CODE_ROOT / "config").is_dir():
        return CODE_ROOT
    return Path.cwd().resolve()


def resource_root() -> Path:
    """Reference databases shared across jobs."""
    value = _env_path(ENV_RESOURCES)
    return value if value is not None else (project_root() / "resources").resolve()


def external_root() -> Path:
    """Out-of-tree material (external datasets, legacy workbenches)."""
    value = _env_path(ENV_EXTERNAL_ROOT)
    return value if value is not None else (project_root() / "external").resolve()


def runtime_root(project: Path | None = None) -> Path:
    """Writable runtime root holding envs, caches, logs and manifests.

    A classic deployment keeps the code in ``<site>/ultimate`` and the runtime
    beside it in ``<site>/ultimate_runtime``, which keeps the code root free of
    machine state. Any other directory name nests the runtime inside it.
    """
    value = _env_path(ENV_RUNTIME_ROOT)
    if value is not None:
        return value
    root = (project or project_root()).resolve()
    if root.name == "ultimate":
        return (root.parent / RUNTIME_DIRNAME).resolve()
    return (root / RUNTIME_DIRNAME).resolve()


def environment_root(project: Path | None = None) -> Path:
    """Managed environments (conda/mamba prefixes) live here."""
    value = _env_path(ENV_ENV_ROOT)
    return value if value is not None else (runtime_root(project) / "envs").resolve()


def cache_root(project: Path | None = None) -> Path:
    """Download and artefact caches shared between jobs."""
    value = _env_path(ENV_CACHE_ROOT)
    return value if value is not None else (runtime_root(project) / "cache").resolve()


def default_platform_root() -> Path:
    """Backward-compatible alias for :func:`project_root`."""
    return project_root()


def default_environment_root(project: Path | None = None) -> Path:
    """Backward-compatible alias for :func:`environment_root`."""
    return environment_root(project)


def default_reference_root() -> Path | None:
    """Optional shared reference tree (genomes, MSigDB, CellChatDB).

    Returns ``None`` when unset so callers can keep their own defaults instead of
    inventing a path that does not exist on this site.
    """
    return _env_path(ENV_REFERENCE_ROOT)
