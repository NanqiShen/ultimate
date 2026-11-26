from __future__ import annotations

from ultimate.paths import project_root

import json
from contextlib import nullcontext
from pathlib import Path

import click

from ultimate.approval_gate import load_production_approval
from ultimate.batch import prepare_batch
from ultimate.batch_status import build_batch_status
from ultimate.config import load_config
from ultimate.constants import PROJECT_TYPES
from ultimate.customer_package import build_customer_package
from ultimate.demo import init_project
from ultimate.delivery_check import run_delivery_check
from ultimate.diagnostics import inspect_environment, inspect_trace
from ultimate.handoff_check import run_handoff_check
from ultimate.job import prepare_job
from ultimate.pipeline import run_pipeline_from_config
from ultimate.plot_style import available_styles, generate_style_review, set_active_style
from ultimate.preflight import run_preflight
from ultimate.production_audit import run_production_audit
from ultimate.raw_upstream import run_raw_upstream_evidence
from ultimate.report import build_report
from ultimate.runtime import RunTrace, apply_thread_budget, thread_budget
from ultimate.reproducibility import export_reproducible_package
from ultimate.singlecell_audit import run_singlecell_audit
from ultimate.triage import run_triage
from ultimate.validation_index import build_validation_index


@click.group()
def main() -> None:
    """Ultimate multi-omics bioinformatics command line interface."""


@main.command("init-project")
@click.option("--type", "project_type", type=click.Choice(PROJECT_TYPES), required=True)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Directory where the project template should be created.",
)
@click.option("--demo-data/--no-demo-data", default=False, show_default=True)
def init_project_command(project_type: str, output_dir: Path, demo_data: bool) -> None:
    manifest = init_project(project_type, output_dir, demo_data=demo_data)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("prepare-job")
@click.option("--config", "config_path", type=click.Path(path_type=Path, exists=True, dir_okay=False), required=True)
@click.option("--job-id", required=True, help="Stable job id under <root>/jobs/<job_id>.")
@click.option("--root", type=click.Path(path_type=Path), default=project_root, show_default=True)
@click.option("--samplesheet", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None)
@click.option("--analysis-request", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None)
@click.option("--run-mode", type=click.Choice(["production", "interactive"]), default="production", show_default=True)
def prepare_job_command(
    config_path: Path,
    job_id: str,
    root: Path,
    samplesheet: Path | None,
    analysis_request: Path | None,
    run_mode: str,
) -> None:
    try:
        manifest = prepare_job(
            config_path=config_path,
            job_id=job_id,
            root=root,
            samplesheet=samplesheet,
            analysis_request=analysis_request,
            run_mode=run_mode,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("prepare-batch")
@click.option("--batch", "batch_path", type=click.Path(path_type=Path, exists=True, dir_okay=False), required=True)
@click.option("--root", type=click.Path(path_type=Path), default=None, help="Override batch root. Defaults to the root declared in the batch file.")
@click.option("--output-dir", type=click.Path(path_type=Path), default=None, help="Where batch scaffold artifacts should be written.")
@click.option("--run-mode", type=click.Choice(["production", "interactive"]), default=None, help="Override run_mode for all jobs unless a job overrides it.")
def prepare_batch_command(batch_path: Path, root: Path | None, output_dir: Path | None, run_mode: str | None) -> None:
    """Scaffold multiple prepared jobs from a batch request without running analysis."""
    try:
        manifest = prepare_batch(batch_path=batch_path, root=root, output_dir=output_dir, run_mode=run_mode)
    except (OSError, ValueError, TypeError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("preflight")
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
)
def preflight_command(config_path: Path) -> None:
    try:
        loaded = load_config(config_path)
        manifest = run_preflight(loaded.raw, write=True)
    except (OSError, TypeError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("doctor")
@click.option("--config", "config_path", type=click.Path(path_type=Path, exists=True, dir_okay=False))
@click.option("--threads", type=click.IntRange(min=1))
@click.option("--module", "module_name", default=None, help="Also check one module's Python/R dependencies and packaged assets.")
def doctor_command(config_path: Path | None, threads: int | None, module_name: str | None) -> None:
    """Check the installation and configuration metadata without running analyses."""
    report = inspect_environment(config_path, threads, module_name)
    click.echo(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "pass":
        raise click.exceptions.Exit(1)


@main.command("debug")
@click.option("--config", "config_path", type=click.Path(path_type=Path, exists=True, dir_okay=False))
@click.option("--trace", "trace_path", type=click.Path(path_type=Path, exists=True, dir_okay=False))
@click.option("--threads", type=click.IntRange(min=1))
def debug_command(config_path: Path | None, trace_path: Path | None, threads: int | None) -> None:
    """Explain the planned CPU budget, or summarize a trace without replaying a run."""
    if not config_path and not trace_path:
        raise click.UsageError("provide --config and/or --trace")
    report = {}
    if config_path:
        report["plan"] = inspect_environment(config_path, threads)
    if trace_path:
        try:
            report["trace"] = inspect_trace(trace_path)
        except (OSError, ValueError) as exc:
            raise click.ClickException(str(exc)) from exc
    click.echo(json.dumps(report, indent=2, ensure_ascii=False))
    if report.get("plan", {}).get("status") == "fail":
        raise click.exceptions.Exit(1)


@main.command("run")
@click.option(
    "--config",
    "config_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
)
@click.option("--production-approval", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None, help="Approved JSON gate file required when the unified run requests production_backend.")
@click.option("--threads", type=click.IntRange(min=1), help="Requested thread count; capped by Slurm and CPU affinity.")
@click.option("--debug-log", type=click.Path(path_type=Path, dir_okay=False), help="New JSONL stage trace; existing files are not overwritten.")
def run_command(config_path: Path, production_approval: Path | None, threads: int | None, debug_log: Path | None) -> None:
    try:
        handle = debug_log.open("x", encoding="utf-8") if debug_log else nullcontext(None)
        with handle as stream:
            manifest = run_pipeline_from_config(config_path, production_approval_path=production_approval,
                                                threads=threads, trace=RunTrace(stream) if stream else None)
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("triage")
@click.option("--request", "request_path", type=click.Path(path_type=Path, exists=True, dir_okay=False), required=True)
@click.option("--output-dir", type=click.Path(path_type=Path), required=True)
def triage_command(request_path: Path, output_dir: Path) -> None:
    """Assess technical readiness without running analysis or creating production artifacts."""
    manifest = run_triage(request_path=request_path, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("report")
@click.option(
    "--run-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    required=True,
)
def report_command(run_dir: Path) -> None:
    manifest = build_report(run_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("export-repro")
@click.option(
    "--run-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    required=True,
    help="Run directory containing run_manifest.json.",
)
@click.option("--checksum-max-mb", type=int, default=256, show_default=True, help="Maximum file size to hash for input checksums.")
def export_repro_command(run_dir: Path, checksum_max_mb: int) -> None:
    manifest = export_reproducible_package(run_dir, checksum_max_bytes=checksum_max_mb * 1024 * 1024)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("delivery-check")
@click.option(
    "--run-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    required=True,
    help="Run directory or prepared job directory to check before delivery.",
)
def delivery_check_command(run_dir: Path) -> None:
    """Run the delivery QA gate for a production run or prepared job."""
    manifest = run_delivery_check(run_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))
    if manifest.get("status") != "ready":
        raise click.ClickException("delivery-check blocked: " + ",".join(manifest.get("blockers") or []))


@main.command("customer-package")
@click.option("--run-dir", type=click.Path(path_type=Path, exists=True, file_okay=False), required=True, help="Production run directory or prepared job directory.")
@click.option("--output-dir", type=click.Path(path_type=Path), default=None, help="Where sanitized customer-facing files should be written.")
def customer_package_command(run_dir: Path, output_dir: Path | None) -> None:
    """Build a sanitized customer-facing package from a production run."""
    manifest = build_customer_package(run_dir=run_dir, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))
    if manifest.get("status") != "ready":
        raise click.ClickException("customer-package blocked: sanitization checks failed")


@main.command("batch-status")
@click.option("--batch-dir", type=click.Path(path_type=Path, exists=True, file_okay=False), required=True, help="Batch output directory or directory containing prepared jobs.")
@click.option("--output-dir", type=click.Path(path_type=Path), default=None, help="Where batch status artifacts should be written. Defaults to --batch-dir.")
@click.option("--job-glob", default=None, help="Optional job directory glob used to limit status scans, for example 'v4_2_*_20260611T000000Z'.")
def batch_status_command(batch_dir: Path, output_dir: Path | None, job_glob: str | None) -> None:
    """Summarize raw-upstream, run, customer-package, and delivery-check status for a batch."""
    manifest = build_batch_status(batch_dir=batch_dir, output_dir=output_dir, job_glob=job_glob)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("handoff-check")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=Path("."),
    show_default=True,
    help="Ultimate project root used for audit output defaults.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where handoff QA artifacts should be written. Defaults to <root>/audits/handoff_latest.",
)
def handoff_check_command(root: Path, output_dir: Path | None) -> None:
    """Check nf-core raw upstream handoff templates without executing upstream workflows."""
    manifest = run_handoff_check(root=root, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))
    if manifest.get("status") != "ready":
        raise click.ClickException("handoff-check blocked: " + ",".join(manifest.get("blockers") or []))


@main.command("check-upstream-acceptance")
@click.option("--package-dir", "package_dirs", type=click.Path(path_type=Path, exists=True, file_okay=False), multiple=True, required=True, help="prepare-upstream package directory (can be specified multiple times).")
@click.option("--downstream-manifest", "downstream_manifests", type=click.Path(path_type=Path, dir_okay=False), multiple=True, required=True, help="run_manifest.json produced by the downstream Ultimate run (can be specified multiple times).")
@click.option("--output", "output_path", type=click.Path(path_type=Path), default=None, help="Where to write the acceptance JSON.")
@click.option("--require-both-routes/--allow-single-route", default=True, show_default=True, help="Require both salmon and starsolo routes for final M6 sign-off.")
def check_upstream_acceptance_command(
    package_dirs: tuple[Path, ...],
    downstream_manifests: tuple[Path, ...],
    output_path: Path | None,
    require_both_routes: bool,
) -> None:
    """Strictly validate an upstream -> import -> run chain (non-zero on any gap)."""
    from ultimate.upstream_acceptance import (
        RouteAcceptance,
        check_downstream_consumed,
        check_expected_files,
        check_matrix_contract,
        check_sample_count_agreement,
        check_salmon_quant_reconciliation,
        check_10x_directory,
        evaluate,
        exit_code,
        write_acceptance,
    )
    import pandas as pd

    # Match each package directory with its corresponding downstream manifest
    packages: list[dict[str, Any]] = []
    for p_dir in package_dirs:
        p_path = p_dir.resolve()
        m_path = p_path / "job_package_manifest.json"
        if not m_path.is_file():
            raise click.ClickException(f"not a prepare-upstream package: {p_path}")
        pkg_data = json.loads(m_path.read_text(encoding="utf-8"))
        packages.append({"dir": p_path, "route": str(pkg_data.get("route") or ""), "data": pkg_data})

    manifest_map: dict[str, Path] = {}
    remaining_manifests = [m.resolve() for m in downstream_manifests]

    # If 1-to-1 matching
    if len(packages) == 1 and len(remaining_manifests) == 1:
        manifest_map[packages[0]["route"]] = remaining_manifests[0]
    else:
        for m_file in remaining_manifests:
            try:
                m_payload = json.loads(m_file.read_text(encoding="utf-8"))
            except Exception:
                continue
            text = json.dumps(m_payload).lower()
            if "salmon" in text or "rnaseq" in text:
                manifest_map["salmon"] = m_file
            elif "starsolo" in text or "scrna" in text or "star" in text:
                manifest_map["starsolo"] = m_file

    route_acceptances: list[RouteAcceptance] = []
    for pkg_item in packages:
        route = pkg_item["route"]
        pkg_dir = pkg_item["dir"]
        artifacts_dir = pkg_dir / "upstream_out"
        acceptance = RouteAcceptance(route=route or "unknown")
        route_acceptances.append(acceptance)

        check_expected_files(acceptance, "package", [
            pkg_dir / "run_upstream.sbatch",
            pkg_dir / "expected_import_project.yaml",
            pkg_dir / "checksum_manifest.tsv",
            pkg_dir / "samplesheet_normalized.tsv",
        ])
        normalized = pkg_dir / "samplesheet_normalized.tsv"

        downstream_manifest = manifest_map.get(route)
        if downstream_manifest is None:
            acceptance.add("downstream_manifest_provided", "fail", f"no downstream manifest matched for route {route}")
            continue

        if route == "salmon":
            matrix = artifacts_dir / "counts_raw.tsv"
            check_expected_files(acceptance, "upstream", [matrix, artifacts_dir / "import_config.yaml"])
            exp_dims = None
            if matrix.is_file() and normalized.is_file():
                expected_samples = check_sample_count_agreement(acceptance, "upstream", matrix, normalized)
                matrix_info = check_matrix_contract(acceptance, "upstream", matrix, expected_samples=expected_samples)
                check_salmon_quant_reconciliation(acceptance, "upstream", matrix, artifacts_dir, expected_samples)
                if matrix_info.get("features") and expected_samples:
                    exp_dims = {"n_features": int(matrix_info["features"]), "n_samples": len(expected_samples)}
            check_downstream_consumed(
                acceptance, "downstream", downstream_manifest,
                expected_inputs=matrix if matrix.is_file() else artifacts_dir,
                expected_dimensions=exp_dims,
            )
        elif route == "starsolo":
            samples = pd.read_csv(normalized, sep="\t")["sample_id"].astype(str).tolist() if normalized.is_file() else []
            expected_sample_dirs: list[Path] = []
            star_dims = None
            for sample in samples:
                s_dir = artifacts_dir / "star" / sample / "filtered"
                expected_sample_dirs.append(s_dir)
                d_info = check_10x_directory(acceptance, f"upstream:{sample}", s_dir)
                if d_info.get("shape") and not star_dims:
                    star_dims = {"n_genes": d_info["shape"][0], "n_cells": d_info["shape"][1]}

            check_downstream_consumed(
                acceptance, "downstream", downstream_manifest,
                expected_inputs=expected_sample_dirs if expected_sample_dirs else artifacts_dir,
                expected_dimensions=star_dims,
            )
        else:
            acceptance.add("route_supported", "fail", f"unsupported route: {route}")

    required_routes = ("salmon", "starsolo") if require_both_routes else ()
    payload = evaluate(route_acceptances, required_routes=required_routes)
    target = (output_path or (packages[0]["dir"] / "M6_ACCEPTANCE.json")).resolve()
    write_acceptance(target, payload)
    click.echo(json.dumps(payload, indent=2, ensure_ascii=False))
    raise click.exceptions.Exit(exit_code(payload))


@main.command("prepare-upstream")
@click.option("--route", type=click.Choice(["salmon", "starsolo"]), required=True, help="Open-source upstream route without nextflow/apptainer.")
@click.option("--samplesheet", type=click.Path(path_type=Path, exists=True, dir_okay=False), required=True, help="FASTQ samplesheet.tsv (sample, fastq_1, fastq_2).")
@click.option("--index", type=click.Path(path_type=Path), default=None, help="salmon index dir (salmon) or STAR genome index dir (starsolo).")
@click.option("--fastq-dir", type=click.Path(path_type=Path, exists=True, file_okay=False), default=None, help="Base directory for relative FASTQ paths.")
@click.option("--output-dir", type=click.Path(path_type=Path), required=True, help="Directory to write the reviewable job package into.")
@click.option("--threads", type=click.IntRange(min=1), help="Requested thread count; capped by Slurm and CPU affinity.")
def prepare_upstream_command(
    route: str,
    samplesheet: Path,
    index: Path | None,
    fastq_dir: Path | None,
    output_dir: Path,
    threads: int | None,
) -> None:
    """Generate a reviewable Slurm upstream job package (never auto-submitted)."""
    from ultimate.prepare_upstream import prepare_upstream_package

    result = prepare_upstream_package(
        route=route,
        samplesheet=samplesheet,
        output_dir=output_dir,
        index=index,
        fastq_dir=fastq_dir,
        threads=threads,
    )
    click.echo(json.dumps(result.to_manifest(), indent=2, ensure_ascii=False))
    if result.status != "ready":
        raise click.ClickException("prepare-upstream blocked: " + (result.blocked_reason or "unknown"))


@main.command("raw-upstream-evidence")
@click.option("--module", "module_name", type=click.Choice(["rnaseq", "scrna"]), required=True)
@click.option("--input-path", type=click.Path(path_type=Path), required=True)
@click.option("--samplesheet", type=click.Path(path_type=Path), default=None)
@click.option("--output-dir", type=click.Path(path_type=Path), required=True)
@click.option("--stage", default=None, help="Optional stage label recorded in the raw upstream manifest.")
@click.option("--tiny-reference", type=click.Path(path_type=Path), default=None, help="Tiny FASTA/reference required for rnaseq_fastq_tiny_counts.")
@click.option("--quant-tool", default=None, help="Quant/count command to require for rnaseq_fastq_tiny_counts: salmon, featureCounts, or subread.")
def raw_upstream_evidence_command(
    module_name: str,
    input_path: Path,
    samplesheet: Path | None,
    output_dir: Path,
    stage: str | None,
    tiny_reference: Path | None,
    quant_tool: str | None,
) -> None:
    """Write lightweight Slurm evidence that raw/semi-raw inputs are importable."""
    manifest = run_raw_upstream_evidence(
        module=module_name,
        input_path=input_path,
        samplesheet=samplesheet,
        output_dir=output_dir,
        stage=stage,
        tiny_reference=tiny_reference,
        quant_tool=quant_tool,
    )
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))
    if manifest.get("status") != "ready":
        raise click.ClickException("raw-upstream-evidence blocked: " + str(manifest.get("blocked_reason") or "unknown"))


@main.command("audit-singlecell")
@click.option(
    "--root",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    default=project_root,
    show_default=True,
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where audit artifacts should be written. Defaults to <root>/audits/singlecell.",
)
def audit_singlecell_command(root: Path, output_dir: Path | None) -> None:
    manifest = run_singlecell_audit(root=root, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("validation-index")
@click.option(
    "--root",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    default=project_root,
    show_default=True,
)
@click.option(
    "--validations-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    default=None,
    help="Directory containing <run>/run_manifest.json files. Defaults to <root>/validations.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where validation index artifacts should be written. Defaults to <root>/reports/validation_index.",
)
def validation_index_command(root: Path, validations_dir: Path | None, output_dir: Path | None) -> None:
    manifest = build_validation_index(root=root, validations_dir=validations_dir, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("audit-production")
@click.option(
    "--root",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    default=project_root,
    show_default=True,
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where production-readiness audit artifacts should be written.",
)
def audit_production_command(root: Path, output_dir: Path | None) -> None:
    manifest = run_production_audit(root=root, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("audit-modules")
@click.option(
    "--root",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    default=None,
    help="Project root to audit. Defaults to the installed ultimate package.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Where module standardization artifacts should be written.",
)
def audit_modules_command(root: Path | None, output_dir: Path) -> None:
    from datetime import datetime, timezone

    import pandas as pd

    from ultimate.module_standardization import build_module_standardization_rows

    output_dir.mkdir(parents=True, exist_ok=True)
    modules_root = root / "src" / "ultimate" / "modules" if root else None
    rows = build_module_standardization_rows(modules_root=modules_root)
    matrix_path = output_dir / "module_standardization_matrix.tsv"
    pd.DataFrame(rows).to_csv(matrix_path, sep="\t", index=False)
    summary = {
        "ready": sum(1 for row in rows if row["overall_status"] == "ready"),
        "partial": sum(1 for row in rows if row["overall_status"] != "ready"),
    }
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "root": str(root.resolve()) if root else "",
        "output_dir": str(output_dir.resolve()),
        "module_count": len(rows),
        "summary": summary,
        "module_standardization_matrix": str(matrix_path.resolve()),
    }
    manifest_path = output_dir / "run_manifest.json"
    manifest["manifest_path"] = str(manifest_path.resolve())
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("audit-tools")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=project_root,
    show_default=True,
    help="Ultimate project root on shared storage.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where tool audit artifacts should be written. Defaults to <root>/audits/tools.",
)
def audit_tools_command(root: Path, output_dir: Path | None) -> None:
    from ultimate.tool_registry import run_audit_tools

    manifest = run_audit_tools(root=root, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("audit-backends")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=project_root,
    show_default=True,
    help="Ultimate project root on shared storage.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where backend audit artifacts should be written. Defaults to <root>/audits/backends.",
)
def audit_backends_command(root: Path, output_dir: Path | None) -> None:
    from ultimate.backend_registry import run_audit_backends

    manifest = run_audit_backends(root=root, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("tool-completeness")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=project_root,
    show_default=True,
    help="Ultimate project root on shared storage.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where V4.1 tool completeness artifacts should be written.",
)
def tool_completeness_command(root: Path, output_dir: Path | None) -> None:
    from ultimate.completeness import run_tool_completeness

    manifest = run_tool_completeness(root=root, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("order-readiness")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=project_root,
    show_default=True,
    help="Ultimate project root on shared storage.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where V4.1 module order-readiness artifacts should be written.",
)
def order_readiness_command(root: Path, output_dir: Path | None) -> None:
    from ultimate.completeness import run_module_order_readiness

    manifest = run_module_order_readiness(root=root, output_dir=output_dir)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("trial-tools")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=project_root,
    show_default=True,
    help="Ultimate project root on shared storage.",
)
@click.option("--batch", required=True)
@click.option("--output-dir", type=click.Path(path_type=Path), default=None)
@click.option("--project-root", type=click.Path(path_type=Path), default=None, help="Directory containing envs/*.yml. Defaults to --root.")
@click.option("--install/--no-install", default=False, show_default=True, help="Run the batch mamba install before smoke checks.")
def trial_tools_command(root: Path, batch: str, output_dir: Path | None, project_root: Path | None, install: bool) -> None:
    from ultimate.tool_registry import available_tool_batches, run_trial_tools

    if batch not in available_tool_batches():
        raise click.BadParameter(f"Unsupported batch {batch!r}; expected one of {available_tool_batches()}")
    manifest = run_trial_tools(root=root, batch=batch, output_dir=output_dir, install=install, project_root=project_root)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("prune-tools")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=project_root,
    show_default=True,
    help="Ultimate project root on shared storage.",
)
@click.option("--output-dir", type=click.Path(path_type=Path), default=None)
@click.option("--yes", is_flag=True, help="Actually run safe cache cleanup commands. Without this, only writes a prune plan.")
def prune_tools_command(root: Path, output_dir: Path | None, yes: bool) -> None:
    from ultimate.tool_registry import run_prune_tools

    manifest = run_prune_tools(root=root, output_dir=output_dir, yes=yes)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("prepare-intake")
@click.option(
    "--root",
    type=click.Path(path_type=Path),
    default=project_root,
    show_default=True,
    help="Ultimate project root on shared storage.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where the customer intake package should be written. Defaults to <root>/intake_packages/latest.",
)
@click.option("--refresh-audit/--no-refresh-audit", default=False, show_default=True)
def prepare_intake_command(root: Path, output_dir: Path | None, refresh_audit: bool) -> None:
    from ultimate.intake import prepare_intake_package

    manifest = prepare_intake_package(root=root, output_dir=output_dir, refresh_audit=refresh_audit)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("styles")
@click.option("--style", "style_key", default="soft_color", show_default=True, help="Style key to render.")
@click.option("--all", "render_all", is_flag=True, help="Render review figures for every registered style.")
@click.option("--output-dir", type=click.Path(path_type=Path), default=None, help="Optional review output directory. Defaults to style_reviews/v3_6/<style>.")
@click.option("--minimal", "minimal", is_flag=True, help="Render sparse review figures with legends/captions/grids off.")
@click.option("--publication", "publication", is_flag=True, help="Render full publication-style review figures. This is the default.")
@click.option("--list", "list_styles", is_flag=True, help="List registered style tokens without rendering review figures.")
def styles_command(style_key: str, render_all: bool, output_dir: Path | None, minimal: bool, publication: bool, list_styles: bool) -> None:
    styles = available_styles()
    options_preset = "minimal" if minimal else "publication"
    if list_styles:
        click.echo(json.dumps(styles, indent=2, ensure_ascii=False))
        return
    if render_all:
        manifests = {}
        output_dir = output_dir or Path("style_reviews") / "v3_6"
        output_dir.mkdir(parents=True, exist_ok=True)
        for key in styles:
            tokens = set_active_style(key)
            manifests[key] = generate_style_review(output_dir / key, style=tokens, options_preset=options_preset)
        click.echo(
            json.dumps(
                {"selected": "all", "available": list(styles), "options_preset": options_preset, "manifests": manifests},
                indent=2,
                ensure_ascii=False,
            )
        )
        return
    output_dir = output_dir or Path("style_reviews") / "v3_6" / style_key
    tokens = set_active_style(style_key)
    manifest = generate_style_review(output_dir, style=tokens, options_preset=options_preset)
    click.echo(json.dumps({"selected": style_key, "available": list(styles), "options_preset": options_preset, **manifest}, indent=2, ensure_ascii=False))


@main.command("create-scrna-demo-inputs")
@click.option("--output-dir", type=click.Path(path_type=Path), required=True)
@click.option("--n-cells", type=int, default=120, show_default=True)
@click.option("--n-genes", type=int, default=90, show_default=True)
@click.option("--seed", type=int, default=17, show_default=True)
def create_scrna_demo_inputs_command(output_dir: Path, n_cells: int, n_genes: int, seed: int) -> None:
    """Create tiny h5ad/10x h5/10x mtx inputs for non-deliverable scRNA MVP checks."""
    from ultimate.scrna_smoke import create_demo_inputs

    manifest = create_demo_inputs(output_dir, n_cells=n_cells, n_genes=n_genes, seed=seed)
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("validate-scrna")
@click.option("--input-path", type=click.Path(path_type=Path, exists=True), required=True)
@click.option("--input-type", type=click.Choice(["h5ad", "10x_h5", "10x_mtx"]), required=True)
@click.option("--output-dir", type=click.Path(path_type=Path), required=True)
@click.option("--samplesheet", type=click.Path(path_type=Path, exists=True), default=None)
@click.option("--max-cells", type=int, default=3000, show_default=True)
@click.option("--random-seed", type=int, default=7, show_default=True)
@click.option("--analysis-level", type=click.Choice(["demo_result", "smoke_backend", "validated_backend", "production_backend"]), default=None)
@click.option("--public-dataset", is_flag=True, help="Mark a real public dataset validation run; never use with generated demo inputs.")
@click.option("--dataset-label", default=None, help="Optional dataset label recorded in the manifest.")
@click.option("--production-approval", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None, help="Approved JSON gate file required for production_backend.")
@click.option("--celltypist-model", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None, help="Optional local CellTypist model. If omitted, CellTypist is skipped without downloading.")
@click.option("--nichenet-resource", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None, help="Optional ligand-target TSV/CSV resource for NicheNet-style communication backend.")
@click.option("--threads", type=click.IntRange(min=1), help="Requested thread count; capped by Slurm and CPU affinity.")
def validate_scrna_command(
    input_path: Path,
    input_type: str,
    output_dir: Path,
    samplesheet: Path | None,
    max_cells: int,
    random_seed: int,
    analysis_level: str | None,
    public_dataset: bool,
    dataset_label: str | None,
    production_approval: Path | None,
    celltypist_model: Path | None,
    nichenet_resource: Path | None,
    threads: int | None,
) -> None:
    """Run the scRNA MVP validation path on h5ad, 10x H5, or 10x MTX input."""
    from ultimate.scrna_smoke import run_scrna_validation

    try:
        approval = load_production_approval(
            production_approval,
            analysis_level=analysis_level,
            input_path=input_path,
            output_dir=output_dir,
        )
        with apply_thread_budget(thread_budget({}, threads)):
            manifest = run_scrna_validation(
                input_path=input_path,
                input_type=input_type,
                output_dir=output_dir,
                samplesheet=samplesheet,
                max_cells=max_cells,
                random_seed=random_seed,
                analysis_level=analysis_level,
                public_dataset=public_dataset,
                dataset_label=dataset_label,
                production_approval=approval,
                celltypist_model=celltypist_model,
                nichenet_resource=nichenet_resource,
            )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    manifest["thread_budget"] = thread_budget({}, threads).as_dict()
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("validate-scrna-pseudobulk-de")
@click.option("--output-dir", type=click.Path(path_type=Path), required=True)
@click.option("--sce-rda", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None)
@click.option("--counts", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None)
@click.option("--design", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None)
@click.option("--rscript", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None)
@click.option("--sce-object", default=None)
@click.option("--assay", default="counts", show_default=True)
@click.option("--donor-column", default="ind", show_default=True)
@click.option("--condition-column", default="stim", show_default=True)
@click.option("--cell-type-column", default="cell", show_default=True)
@click.option("--cell-type", required=True)
@click.option("--control", required=True)
@click.option("--case", "case_", required=True)
@click.option("--multiplet-column", default=None)
@click.option("--multiplet-keep", default=None)
@click.option("--min-cells", type=int, default=1, show_default=True)
@click.option("--min-pairs", type=int, default=4, show_default=True)
@click.option("--analysis-level", type=click.Choice(["demo_result", "smoke_backend", "validated_backend"]), default="validated_backend", show_default=True)
@click.option("--robust/--no-robust", default=True, show_default=True)
@click.option("--threads", type=click.IntRange(min=1), help="Requested thread count; capped by Slurm and CPU affinity.")
def validate_scrna_pseudobulk_de_command(
    output_dir: Path,
    sce_rda: Path | None,
    counts: Path | None,
    design: Path | None,
    rscript: Path | None,
    sce_object: str | None,
    assay: str,
    donor_column: str,
    condition_column: str,
    cell_type_column: str,
    cell_type: str,
    control: str,
    case_: str,
    multiplet_column: str | None,
    multiplet_keep: str | None,
    min_cells: int,
    min_pairs: int,
    analysis_level: str,
    robust: bool,
    threads: int | None,
) -> None:
    """Run paired edgeR pseudobulk DE from SCE or preaggregated pseudobulk inputs."""
    from ultimate.scrna_pseudobulk_de import run_scrna_pseudobulk_de

    budget = thread_budget({}, threads)
    try:
        manifest = run_scrna_pseudobulk_de(
            output_dir=output_dir,
            sce_rda=sce_rda,
            counts=counts,
            design=design,
            rscript=rscript,
            sce_object=sce_object,
            assay=assay,
            donor_column=donor_column,
            condition_column=condition_column,
            cell_type_column=cell_type_column,
            cell_type=cell_type,
            control=control,
            case=case_,
            multiplet_column=multiplet_column,
            multiplet_keep=multiplet_keep,
            min_cells=min_cells,
            min_pairs=min_pairs,
            analysis_level=analysis_level,
            robust=robust,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        raise click.ClickException(str(exc)) from exc
    manifest["thread_budget"] = budget.as_dict()
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


@main.command("annotate-scrna-celltypist")
@click.option("--input-h5ad", type=click.Path(path_type=Path, exists=True, dir_okay=False), required=True)
@click.option("--output-dir", type=click.Path(path_type=Path), required=True)
@click.option("--model", "model_path", type=click.Path(path_type=Path, exists=True, dir_okay=False), required=True)
@click.option("--analysis-level", type=click.Choice(["demo_result", "smoke_backend", "validated_backend", "production_backend"]), default="validated_backend", show_default=True)
@click.option("--cluster-key", default="cluster", show_default=True)
@click.option("--embedding-h5ad", type=click.Path(path_type=Path, exists=True, dir_okay=False), default=None, help="Optional h5ad with UMAP/cluster metadata aligned by cell id.")
@click.option("--majority-voting/--no-majority-voting", default=True, show_default=True)
@click.option("--uncertainty-threshold", type=float, default=0.5, show_default=True, help="Exploratory max-probability threshold used only to flag uncertain cells.")
@click.option("--threads", type=click.IntRange(min=1), help="Requested thread count; capped by Slurm and CPU affinity.")
def annotate_scrna_celltypist_command(
    input_h5ad: Path,
    output_dir: Path,
    model_path: Path,
    analysis_level: str,
    cluster_key: str,
    embedding_h5ad: Path | None,
    majority_voting: bool,
    uncertainty_threshold: float,
    threads: int | None,
) -> None:
    """Run CellTypist annotation on an existing full-gene h5ad without rerunning the scRNA MVP pipeline."""
    from ultimate.scrna_smoke import run_celltypist_annotation

    budget = thread_budget({}, threads)
    try:
        manifest = run_celltypist_annotation(
            input_h5ad=input_h5ad,
            output_dir=output_dir,
            model_path=model_path,
            analysis_level=analysis_level,
            cluster_key=cluster_key,
            embedding_h5ad=embedding_h5ad,
            majority_voting=majority_voting,
            uncertainty_threshold=uncertainty_threshold,
        )
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise click.ClickException(str(exc)) from exc
    manifest["thread_budget"] = budget.as_dict()
    click.echo(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
