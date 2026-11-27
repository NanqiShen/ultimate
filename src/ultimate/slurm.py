"""Cluster submission entry points.

Sites differ in how jobs reach the scheduler: some clusters ship a thin
``sbatch`` wrapper, others expose the native binary only. Ultimate always emits
the standard ``sbatch`` invocation and lets a site override the binary through
``ULTIMATE_SBATCH`` instead of hard-coding a wrapper name.
"""

from __future__ import annotations

import os

ENV_SBATCH = "ULTIMATE_SBATCH"
DEFAULT_SBATCH = "sbatch"


def sbatch_command() -> str:
    """The sbatch binary this site uses, e.g. ``sbatch`` or an absolute path."""
    return os.environ.get(ENV_SBATCH) or DEFAULT_SBATCH


def sbatch_line(script: str, *, args: str = "") -> str:
    """Build a copy-pasteable submission command for ``script``."""
    command = f"{sbatch_command()} {script}"
    return f"{command} {args}".rstrip() if args else command
