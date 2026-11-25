"""Metadata-only pre-run checks and inspection of explicitly supplied traces."""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys

import yaml
from pathlib import Path

from ultimate import __version__
from ultimate.config import enabled_modules, load_config
from ultimate.paths import project_root, resource_root
from ultimate.runtime import thread_budget
from ultimate.template_resources import find_template_dir

_CORE = ("click", "jinja2", "matplotlib", "numpy", "pandas", "yaml", "seaborn", "threadpoolctl")
_PHASES = ("config", "run", "approval", "preflight", "raw_qc", "module", "finalize")


def inspect_environment(config_path: Path | None = None, threads: int | None = None) -> dict:
    checks = []
    config = {}
    if config_path is not None:
        try:
            config = load_config(config_path).raw
        except (OSError, ValueError, TypeError, AttributeError, yaml.YAMLError) as exc:
            checks.append({"name": "config", "status": "fail", "detail": type(exc).__name__})
        else:
            checks.append({"name": "config", "status": "pass"})
            target = Path(config["project"]["output_dir"])
            ancestor = target
            while not ancestor.exists() and not ancestor.is_symlink() and ancestor != ancestor.parent:
                ancestor = ancestor.parent
            writable = ancestor.is_dir() and os.access(ancestor, os.W_OK | os.X_OK)
            checks.append({"name": "output_location", "status": "pass" if writable else "fail",
                           "path": str(target), "checked_ancestor": str(ancestor)})
            sheet = config.get("samples") or {}
            if isinstance(sheet, dict) and sheet.get("samplesheet"):
                p = Path(sheet["samplesheet"])
                checks.append({"name": "samplesheet_metadata", "status": "pass" if p.is_file() and os.access(p, os.R_OK) else "fail"})
    checks.append({"name": "python", "status": "pass" if sys.version_info >= (3, 11) else "fail", "version": sys.version.split()[0]})
    for name in _CORE:
        checks.append({"name": f"dependency:{name}", "status": "pass" if importlib.util.find_spec(name) else "fail"})
    for name in ("intake", "samples"):
        checks.append({"name": f"templates:{name}", "status": "pass" if find_template_dir(name) else "fail"})
    for tool in ("sbatch", "Rscript"):
        checks.append({"name": f"optional:{tool}", "status": "pass" if shutil.which(tool) else "warn"})
    try:
        budget = thread_budget(config, threads).as_dict()
    except ValueError as exc:
        checks.append({"name": "threads", "status": "fail", "detail": str(exc)})
        budget = None
    return {"status": "fail" if any(c["status"] == "fail" for c in checks) else "pass",
            "version": __version__, "project_root": str((config.get("project") or {}).get("server_root") or project_root()),
            "resource_root": str(resource_root()), "threads": budget, "modules": enabled_modules(config),
            "checks": checks, "scope": "metadata_only_not_production_approval"}


def inspect_trace(path: Path) -> dict:
    """Summarize known stages without echoing arbitrary log text."""
    counts = {}
    last_stage = None
    failed = None
    limit = 2 * 1024 * 1024
    with path.open("rb") as f:
        payload = f.read(limit + 1)
    if len(payload) > limit:
        raise ValueError("trace exceeds 2 MiB")
    for number, raw in enumerate(payload.splitlines(), 1):
        try:
            event = json.loads(raw)
        except (UnicodeError, ValueError) as exc:
            raise ValueError(f"invalid JSONL event at line {number}") from exc
        if not isinstance(event, dict) or event.get("stage") not in _PHASES or event.get("status") not in ("started", "completed", "failed"):
            raise ValueError(f"unknown trace event at line {number}")
        last_stage = event["stage"]
        status = event["status"]
        counts[status] = counts.get(status, 0) + 1
        if status == "failed" and failed is None:
            failed = last_stage
    return {"trace_events": counts, "first_failed_stage": failed, "last_stage": last_stage}
