"""Shared path defaults for tools executed directly from their checkout."""
import os
from pathlib import Path


def project_root():
    return Path(os.environ.get("ULTIMATE_ROOT") or Path(__file__).resolve().parents[1]).expanduser().resolve()


def external_root():
    return Path(os.environ.get("ULTIMATE_EXTERNAL_ROOT") or project_root() / "external").expanduser().resolve()


def external_root_0518():
    return external_root() / "0518"
