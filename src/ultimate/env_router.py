"""Cross-environment routing for interpreters and packaged assets.

XieLab keeps every runtime artefact under ``ULTIMATE_RUNTIME_ROOT`` and does not
export conda, mamba, micromamba or Rscript on the default ``$PATH``. Code that
calls ``shutil.which("Rscript")`` therefore cannot see the isolated R 4.4
environment that actually owns ``limma``/``DESeq2``/``edgeR``.

This module is the single place that answers two questions:

* which ``Rscript`` should run a given module's backend, and
* where does a file shipped with the project live (source checkout, ``ULTIMATE_ROOT``
  or ``sys.prefix/share/ultimate`` after a wheel install).

Resolution is always recorded, never guessed: callers receive ``None`` (R) or a
raised ``FileNotFoundError`` (assets) instead of a plausible-looking wrong path.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

from ultimate.paths import environment_root, project_root

ENV_RSCRIPT = "ULTIMATE_RSCRIPT"
ENV_R_ENV_PREFIX = "ULTIMATE_R_ENV_"

DEFAULT_R_ENVS: tuple[str, ...] = ("ultimate-bio-r", "ultimate-rnaseq", "ultimate-scrna-r")


def module_env_variable(module_name: str) -> str:
    """Environment variable that pins the R environment for one module."""
    normalised = "".join(char if char.isalnum() else "_" for char in str(module_name)).upper()
    return f"{ENV_R_ENV_PREFIX}{normalised}"


def _usable_rscript(path: Path | str | None) -> Path | None:
    if not path:
        return None
    candidate = Path(str(path)).expanduser()
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate.resolve()
    return None


def _configured_rscript_values(module_cfg: dict | None) -> list[str]:
    if not isinstance(module_cfg, dict):
        return []
    candidates: list[object] = []
    de_backend = module_cfg.get("de_backend")
    if isinstance(de_backend, dict):
        candidates.extend([de_backend.get("rscript"), de_backend.get("rscript_path")])
    candidates.extend([module_cfg.get("rscript"), module_cfg.get("rscript_path")])
    return [str(value) for value in candidates if value]


def _configured_rscript(module_cfg: dict | None) -> Path | None:
    for value in _configured_rscript_values(module_cfg):
        found = _usable_rscript(value)
        if found is not None:
            return found
    return None


def _env_rscript(module_name: str) -> Path | None:
    for name in (module_env_variable(module_name), ENV_RSCRIPT):
        found = _usable_rscript(os.environ.get(name))
        if found is not None:
            return found
    return None


def resolve_rscript(
    module_name: str,
    module_cfg: dict | None = None,
    preferred_envs: tuple[str, ...] = DEFAULT_R_ENVS,
) -> Path | None:
    """Locate the ``Rscript`` for ``module_name``.

    Priority, highest first:

    1. explicit configuration (``module_cfg["rscript"]`` or
       ``module_cfg["de_backend"]["rscript"]``),
    2. ``ULTIMATE_R_ENV_<MODULE>`` then the global ``ULTIMATE_RSCRIPT``,
    3. ``environment_root()/<env>/bin/Rscript`` for each preferred environment,
    4. ``Rscript`` on ``$PATH``.

    Returns ``None`` when nothing usable exists; callers must record a skipped
    backend instead of silently running without R.
    """
    configured = _configured_rscript(module_cfg)
    if configured is not None:
        return configured
    declared = _configured_rscript_values(module_cfg)
    if declared:
        # An explicit but unusable path is a configuration error; silently
        # substituting another R would run the wrong interpreter.
        raise ValueError(
            f"configured Rscript is not an executable file: {declared[0]}"
        )
    from_env = _env_rscript(module_name)
    if from_env is not None:
        return from_env
    for env_name in preferred_envs:
        found = _usable_rscript(environment_root() / env_name / "bin" / "Rscript")
        if found is not None:
            return found
    found = shutil.which("Rscript")
    return _usable_rscript(found) if found else None


def describe_rscript(
    module_name: str,
    module_cfg: dict | None = None,
    preferred_envs: tuple[str, ...] = DEFAULT_R_ENVS,
) -> dict[str, object]:
    """Explain a :func:`resolve_rscript` result for reports and ``doctor``."""
    configured = _configured_rscript(module_cfg)
    if configured is not None:
        return {"status": "ready", "path": str(configured), "source": "module_config"}
    declared = _configured_rscript_values(module_cfg)
    if declared:
        return {"status": "missing", "path": "", "source": "module_config_invalid",
                "searched": declared[0]}
    from_env = _env_rscript(module_name)
    if from_env is not None:
        source = module_env_variable(module_name) if os.environ.get(module_env_variable(module_name)) else ENV_RSCRIPT
        return {"status": "ready", "path": str(from_env), "source": source}
    env_root = environment_root()
    for env_name in preferred_envs:
        candidate = env_root / env_name / "bin" / "Rscript"
        if _usable_rscript(candidate) is not None:
            return {"status": "ready", "path": str(candidate.resolve()), "source": f"environment_root:{env_name}"}
    on_path = shutil.which("Rscript")
    if on_path:
        return {"status": "ready", "path": str(Path(on_path).resolve()), "source": "PATH"}
    return {
        "status": "missing",
        "path": "",
        "source": "not_found",
        "searched": ";".join(
            [
                f"module_cfg.rscript",
                module_env_variable(module_name),
                ENV_RSCRIPT,
                *[str(env_root / env / "bin" / "Rscript") for env in preferred_envs],
                "PATH:Rscript",
            ]
        ),
    }


def packaged_asset_candidates(*parts: str) -> list[Path]:
    """Asset locations in lookup order for source checkouts and wheel installs."""
    suffix = Path(*parts)
    candidates: list[Path] = []
    override = os.environ.get("ULTIMATE_ASSET_ROOT")
    if override:
        candidates.append(Path(override).expanduser() / suffix)
    root_override = os.environ.get("ULTIMATE_ROOT")
    if root_override:
        candidates.append(Path(root_override).expanduser() / suffix)
    code_root = Path(__file__).resolve().parents[2]
    candidates.append(code_root / suffix)
    candidates.append(Path(sys.prefix) / "share" / "ultimate" / suffix)
    return candidates


def resolve_packaged_asset(*parts: str) -> Path:
    """Find a file that ships with Ultimate.

    Lookup order honours explicit site overrides first, then the source checkout,
    then the ``share/ultimate`` tree produced by a wheel install.
    """
    searched = packaged_asset_candidates(*parts)
    for candidate in searched:
        if candidate.is_file():
            return candidate.resolve()
    joined = f"{'/'.join(parts)}"
    listing = ";".join(str(path) for path in searched)
    raise FileNotFoundError(f"Packaged asset not found: {joined}. Searched: {listing}")


def find_packaged_asset(*parts: str) -> Path | None:
    """Non-raising variant of :func:`resolve_packaged_asset`."""
    for candidate in packaged_asset_candidates(*parts):
        if candidate.is_file():
            return candidate.resolve()
    return None


def packaged_asset_status(*parts: str) -> dict[str, str]:
    found = find_packaged_asset(*parts)
    return {
        "status": "ready" if found else "missing",
        "path": str(found or ""),
        "searched": ";".join(str(candidate) for candidate in packaged_asset_candidates(*parts)),
    }


def scripts_r_dir() -> Path | None:
    """Directory holding the bundled ``*.R`` backends."""
    for candidate in packaged_asset_candidates("scripts", "R"):
        if candidate.is_dir():
            return candidate.resolve()
    return None


def default_config_path() -> Path | None:
    return find_packaged_asset("config", "default.yaml")


def ultimate_run_sbatch_path() -> Path | None:
    return find_packaged_asset("slurm", "ultimate_run.sbatch")


__all__ = [
    "DEFAULT_R_ENVS",
    "ENV_RSCRIPT",
    "ENV_R_ENV_PREFIX",
    "default_config_path",
    "describe_rscript",
    "find_packaged_asset",
    "module_env_variable",
    "packaged_asset_candidates",
    "packaged_asset_status",
    "resolve_packaged_asset",
    "resolve_rscript",
    "scripts_r_dir",
    "ultimate_run_sbatch_path",
]
