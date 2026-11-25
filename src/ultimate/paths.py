"""Deployment paths, resolved at runtime rather than pinned to a cluster."""

from __future__ import annotations

import os
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[2]


def project_root() -> Path:
    """Execution root. Supports source checkouts and explicit site overrides."""
    value = os.environ.get("ULTIMATE_ROOT")
    if value:
        return Path(value).expanduser().resolve()
    # Installed wheels live under site-packages, not a project checkout.
    if (CODE_ROOT / "pyproject.toml").is_file() and (CODE_ROOT / "config").is_dir():
        return CODE_ROOT
    return Path.cwd().resolve()


def resource_root() -> Path:
    value = os.environ.get("ULTIMATE_RESOURCES")
    return Path(value).expanduser().resolve() if value else project_root() / "resources"


def external_root() -> Path:
    value = os.environ.get("ULTIMATE_EXTERNAL_ROOT")
    return Path(value).expanduser().resolve() if value else project_root() / "external"
