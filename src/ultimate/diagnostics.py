"""Metadata-only pre-run checks and inspection of explicitly supplied traces."""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys

import yaml
from pathlib import Path

from ultimate import __version__
from ultimate.bulk import BULK_PYTHON_REQUIREMENTS
from ultimate.config import enabled_modules, load_config
from ultimate.constants import MODULE_SPECS
from ultimate.env_router import (
    default_config_path,
    describe_rscript,
    resolve_rscript,
    scripts_r_dir,
    ultimate_run_sbatch_path,
)
from ultimate.paths import project_root, resource_root
from ultimate.runtime import thread_budget
from ultimate.template_resources import find_template_dir

_CORE = ("click", "jinja2", "matplotlib", "numpy", "pandas", "yaml", "seaborn", "threadpoolctl")
_PHASES = ("config", "run", "approval", "preflight", "raw_qc", "module", "finalize")


def inspect_environment(config_path: Path | None = None, threads: int | None = None, module_name: str | None = None) -> dict:
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
    report = {"status": "fail" if any(c["status"] == "fail" for c in checks) else "pass",
              "version": __version__, "project_root": str((config.get("project") or {}).get("server_root") or project_root()),
              "resource_root": str(resource_root()), "threads": budget, "modules": enabled_modules(config),
              "checks": checks, "scope": "metadata_only_not_production_approval"}
    if module_name:
        module_report = inspect_module(module_name, config or None)
        report["module"] = module_report
        if module_report["status"] == "fail":
            report["status"] = "fail"
    return report


def inspect_module(module_name: str, config: dict | None = None) -> dict:
    """Check one module's optional Python/R dependencies and packaged assets.

    Nothing here runs an analysis or an R process; absent pieces are reported so
    a run can be blocked *before* it silently degrades to a fallback.
    """
    spec = MODULE_SPECS.get(module_name)
    if spec is None:
        return {"name": module_name, "status": "fail", "reason": "unknown_module",
                "checks": [{"name": f"module:{module_name}", "status": "fail", "detail": "unknown module"}],
                "known_modules": sorted(MODULE_SPECS)}

    module_cfg: dict = {}
    if isinstance(config, dict):
        modules = config.get("modules") or {}
        if isinstance(modules, dict) and isinstance(modules.get(module_name), dict):
            module_cfg = modules[module_name]

    checks: list[dict] = []
    python_requirements = BULK_PYTHON_REQUIREMENTS.get(module_name, ())
    for package in sorted(set(python_requirements)):
        available = importlib.util.find_spec(package) is not None
        checks.append({"name": f"python:{package}", "status": "pass" if available else "warn",
                       "detail": "" if available else "optional_python_package_missing"})

    rscript = describe_rscript(module_name, module_cfg or None)
    checks.append({"name": "rscript", "status": "pass" if rscript["status"] == "ready" else "warn",
                   "path": str(rscript.get("path") or ""), "source": str(rscript.get("source") or "")})

    scripts_dir = scripts_r_dir()
    checks.append({"name": "scripts:R", "status": "pass" if scripts_dir else "fail",
                   "path": str(scripts_dir or "")})

    r_packages = _available_r_packages(resolve_rscript(module_name, module_cfg or None), spec.optional_r_packages)
    if r_packages is None:
        checks.append({"name": "r_packages", "status": "warn", "detail": "rscript_unavailable",
                       "declared": list(spec.optional_r_packages)})
    else:
        for package, available in r_packages.items():
            checks.append({"name": f"r:{package}", "status": "pass" if available else "warn",
                           "detail": "" if available else "optional_r_package_missing"})

    wrapper = ultimate_run_sbatch_path()
    checks.append({"name": "slurm:ultimate_run.sbatch", "status": "pass" if wrapper else "fail",
                   "path": str(wrapper or "")})
    config_file = default_config_path()
    checks.append({"name": "config:default.yaml", "status": "pass" if config_file else "warn",
                   "path": str(config_file or "")})

    return {
        "name": module_name,
        "title_cn": spec.title_cn,
        "status": "fail" if any(check["status"] == "fail" for check in checks) else "pass",
        "checks": checks,
        "rscript": rscript,
        "notes": ["metadata_only: no analysis or R process was executed"],
    }


def _available_r_packages(rscript: Path | None, packages: tuple[str, ...]) -> dict[str, bool] | None:
    if rscript is None or not packages:
        return None
    quoted = ",".join(f'"{package}"' for package in packages)
    expression = (
        f"pkgs <- c({quoted}); "
        'for (p in pkgs) { cat(p, requireNamespace(p, quietly=TRUE), sep="\\t"); cat("\\n") }'
    )
    try:
        completed = subprocess.run([str(rscript), "-e", expression], check=False, text=True, capture_output=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    found: dict[str, bool] = {package: False for package in packages}
    for line in completed.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 2 and parts[0] in found:
            found[parts[0]] = parts[1].upper() == "TRUE"
    return found


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
