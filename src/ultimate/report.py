from __future__ import annotations

import base64
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape


TEMPLATE_DIR = Path(__file__).resolve().parent / "templates"


def build_report(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing run manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    reports_dir = run_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    module_visuals = _load_module_visuals(manifest, run_dir)
    context = {
        "manifest": manifest,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "modules": manifest.get("modules", []),
        "preflight": manifest.get("preflight", {}),
        "delivery_gate": manifest.get("delivery_gate", {}),
        "module_qc_manifests": _load_module_qc_manifests(manifest),
        "module_visuals": module_visuals,
        "summary": _build_summary(manifest, module_visuals),
        "embed_images": bool((manifest.get("report") or {}).get("embed_images", False)),
    }
    html = _render("report.html.j2", context)
    methods = _render("methods.md.j2", context)
    html_path = reports_dir / "report.html"
    methods_path = reports_dir / "methods.md"
    report_manifest_path = reports_dir / "report_manifest.json"
    html_path.write_text(html, encoding="utf-8")
    methods_path.write_text(methods, encoding="utf-8")
    report_manifest = {
        "generated_at": context["generated_at"],
        "html_report": str(html_path),
        "methods_md": str(methods_path),
        "source_manifest": str(manifest_path),
        "delivery_gate": manifest.get("delivery_gate", {}),
    }
    report_manifest_path.write_text(json.dumps(report_manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    report_manifest["manifest_path"] = str(report_manifest_path)
    return report_manifest


def _render(template_name: str, context: dict[str, Any]) -> str:
    environment = Environment(
        loader=FileSystemLoader(str(TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=("html", "xml")),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    return environment.get_template(template_name).render(**context)


def _load_module_qc_manifests(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    loaded: dict[str, dict[str, Any]] = {}
    for module in manifest.get("modules", []):
        if not isinstance(module, dict):
            continue
        path_value = ((module.get("artifacts") or {}).get("tables") or {}).get("module_qc_manifest")
        if not path_value:
            continue
        path = Path(str(path_value))
        if not path.exists():
            continue
        try:
            loaded[str(path)] = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            loaded[str(path)] = {"status": "invalid_json"}
    return loaded



# Tables whose first 10 rows are shown as a preview in the HTML report.
TOP_TABLE_KEYS = (
    "differential_results", "de_results", "limma_de_results", "dmp_limma_results",
    "qc_metrics", "survival_results", "pseudobulk_de_results", "marker_genes",
    "enrichment_results", "clinical_feature_associations",
)


def _build_summary(manifest: dict[str, Any], visuals: dict[str, dict[str, Any]]) -> dict[str, Any]:
    modules = manifest.get("modules") or []
    gate = manifest.get("delivery_gate") or {}
    warnings: list[str] = []
    for module in modules:
        if not isinstance(module, dict):
            continue
        level = str(module.get("analysis_level", ""))
        if level != "production_backend":
            warnings.append(f"{module.get('module', 'module')}: analysis_level={level or 'not_recorded'}")
    for module_visual in visuals.values():
        for note in module_visual.get("warnings", []):
            warnings.append(note)
    threads = manifest.get("thread_budget") or {}
    return {
        "n_modules": len(modules),
        "n_samples": _count_samples(modules),
        "gate_status": gate.get("status", "not_recorded"),
        "delivery_allowed": bool(gate.get("delivery_allowed", False)),
        "threads": threads.get("threads", "NA"),
        "warnings": warnings[:20],
    }


def _count_samples(modules: list[Any]) -> int:
    values: set[str] = set()
    for module in modules:
        if not isinstance(module, dict):
            continue
        n = module.get("n_samples")
        if isinstance(n, int):
            values.add(str(n))
    return max((int(v) for v in values), default=0)


def _relative_to_run(path_text: str, run_dir: Path) -> str:
    try:
        path = Path(str(path_text)).resolve()
        return str(path.relative_to(run_dir.resolve())).replace("\\", "/")
    except (ValueError, OSError):
        return str(path_text)


def _read_top_table(path: Path) -> list[dict[str, str]]:
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str).head(10)
    except Exception:
        return []
    return frame.fillna("").to_dict(orient="records")


def _data_uri(path: Path) -> str:
    mime = "image/png" if path.suffix.lower() == ".png" else "image/svg+xml" if path.suffix.lower() == ".svg" else "application/octet-stream"
    try:
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    except OSError:
        return ""
    return f"data:{mime};base64,{encoded}"


def _load_module_visuals(manifest: dict[str, Any], run_dir: Path) -> dict[str, dict[str, Any]]:
    """Collect real figures and Top-10 tables for every module."""
    result: dict[str, dict[str, Any]] = {}
    for module in manifest.get("modules", []):
        if not isinstance(module, dict):
            continue
        name = str(module.get("module", "unknown"))
        artifacts = module.get("artifacts") or {}
        figures: list[dict[str, str]] = []
        for key, value in (artifacts.get("figures") or {}).items():
            if not value:
                continue
            path = Path(str(value))
            if not path.exists():
                continue
            figures.append({
                "key": str(key),
                "src": _relative_to_run(value, run_dir),
                "data_uri": _data_uri(path),
                "label": str(key).replace("_", " "),
            })
        tables: list[dict[str, Any]] = []
        used_paths: set[str] = set()
        for key, value in (artifacts.get("tables") or {}).items():
            if key not in TOP_TABLE_KEYS or not value or str(value) in used_paths:
                continue
            path = Path(str(value))
            if not path.exists():
                continue
            used_paths.add(str(value))
            rows = _read_top_table(path)
            if rows:
                tables.append({"key": str(key), "columns": list(rows[0].keys()), "rows": rows})
        status_pieces = []
        for key in ("wgcna_status", "survival_status"):
            value = module.get(key)
            if value and value != "ready":
                status_pieces.append(f"{key}={value}")
        result[name] = {
            "figures": figures,
            "tables": tables,
            "warnings": status_pieces,
        }
    return result
