"""Generate a reviewable Slurm job package for open-source upstream quantification.

XieLab's compute_pro has no nextflow/apptainer/cellranger, but micromamba can build
an isolated upstream environment (fastp/salmon/STAR/...). This module turns a FASTQ
samplesheet plus a reference index into a self-contained package an operator can
inspect and submit manually:

    run_upstream.sbatch        # the job script (never auto-submitted)
    expected_import_project.yaml
    checksum_manifest.tsv
    job_package_manifest.json

Guarantees: raw FASTQ is only ever read (never copied or modified), threads are
bound to ``$SLURM_CPUS_PER_TASK``/``thread_budget``, and a missing samplesheet,
missing FASTQ or missing reference index blocks the package with a concrete reason.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ultimate.env_router import resolve_packaged_asset
from ultimate.runtime import thread_budget

ROUTES: dict[str, dict[str, str]] = {
    "salmon": {
        "module": "rnaseq",
        "template": "upstream_bulk_rnaseq_salmon.sbatch",
        "requires": "salmon_index",
        "output_hint": "counts_raw.tsv",
    },
    "starsolo": {
        "module": "scrna",
        "template": "upstream_scrna_starsolo.sbatch",
        "requires": "star_index",
        "output_hint": "filtered/",
    },
}

FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

# Index validation markers per route (a directory without these is not usable).
INDEX_MARKERS: dict[str, tuple[str, ...]] = {
    "salmon": ("info.json", "versionInfo.json"),
    "starsolo": ("SAindex", "Genome"),
}


@dataclass
class PrepareUpstreamResult:
    status: str
    route: str
    output_dir: Path
    blocked_reason: str = ""
    blockers: list[str] = field(default_factory=list)
    artifacts: dict[str, str] = field(default_factory=dict)
    n_samples: int = 0
    threads: int = 0

    def to_manifest(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "route": self.route,
            "blocked_reason": self.blocked_reason,
            "blockers": self.blockers,
            "output_dir": str(self.output_dir),
            "n_samples": self.n_samples,
            "threads": self.threads,
            "artifacts": self.artifacts,
        }


def _read_samplesheet(path: Path) -> tuple[list[dict[str, str]], str]:
    """Parse the FASTQ samplesheet; return (rows, error)."""
    if not path.is_file():
        return [], f"missing_samplesheet:{path}"
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        return [], f"unreadable_samplesheet:{type(exc).__name__}"
    if not rows:
        return [], "empty_samplesheet:no_data_rows"
    normalized: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=2):
        clean = {str(key).strip(): str(value or "").strip() for key, value in row.items() if key}
        sample = clean.get("sample") or clean.get("sample_id")
        fq1 = (
            clean.get("fastq_1")
            or clean.get("fq1")
            or clean.get("r1")
            or clean.get("fastq_r1")
            or clean.get("read1")
        )
        if not sample:
            return [], f"samplesheet_missing_sample_at_row:{index}"
        if not fq1:
            return [], f"samplesheet_missing_fastq_1_at_row:{index}"
        clean["sample_id"] = sample
        clean["fastq_1"] = fq1
        normalized.append(clean)
    return normalized, ""


def _resolve_fastq(value: str, base_dir: Path) -> Path:
    candidate = Path(value).expanduser()
    return candidate if candidate.is_absolute() else (base_dir / candidate).resolve()


def _is_fastq(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in FASTQ_SUFFIXES)


def _validate_inputs(
    rows: list[dict[str, str]],
    fastq_dir: Path,
    index: Path,
    route: str,
) -> list[str]:
    """Return blocker strings for read-only input validation."""
    blockers: list[str] = []
    if not index.is_dir():
        blockers.append(f"missing_index_directory:{index}")
    else:
        markers = INDEX_MARKERS.get(route, ())
        if markers and not any((index / marker).exists() for marker in markers):
            blockers.append(f"index_not_usable_for_{route}:missing_any_of_{'|'.join(markers)}")
    for row in rows:
        for key in ("fastq_1", "fastq_2"):
            value = row.get(key, "")
            if not value or value == "-":
                continue
            path = _resolve_fastq(value, fastq_dir)
            if not path.is_file():
                blockers.append(f"missing_fastq:{row['sample_id']}:{key}:{path}")
                continue
            if not os.access(path, os.R_OK):
                blockers.append(f"unreadable_fastq:{row['sample_id']}:{key}")
                continue
            if not _is_fastq(path):
                blockers.append(f"not_a_fastq_file:{row['sample_id']}:{key}:{path.name}")
    # Raw data is strictly read-only: refuse an output dir inside the raw tree.
    return blockers


def _is_inside(candidate: Path, parent: Path) -> bool:
    """True when ``candidate`` is ``parent`` itself or nested below it."""
    try:
        candidate_resolved = candidate.resolve()
        parent_resolved = parent.resolve()
    except OSError:
        return False
    return candidate_resolved == parent_resolved or parent_resolved in candidate_resolved.parents


def _checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _render_run_upstream(
    *,
    route: str,
    samplesheet: Path,
    fastq_dir: Path,
    index: Path,
    out_dir: Path,
    threads: int,
    template_body: str,
    template_name: str,
) -> str:
    body = template_body.split("\n")
    sbatch_lines = [line for line in body if line.startswith("#SBATCH")]
    # Keep the requested CPU count consistent with the thread budget so the job
    # never over- or under-subscribes relative to what the pipeline assumes.
    rebound: list[str] = []
    for line in sbatch_lines:
        if line.startswith("#SBATCH --cpus-per-task="):
            rebound.append(f"#SBATCH --cpus-per-task={threads}")
        else:
            rebound.append(line)
    sbatch_lines = rebound
    rest = [line for line in body if not line.startswith("#SBATCH") and not line.startswith("#!")]
    # Slurm only honours #SBATCH directives inside the leading comment block, so
    # the shebang + directives must precede every executable statement.
    preamble = [
        "#!/usr/bin/env bash",
        f"# Generated by `ultimate prepare-upstream` -- route={route}",
        f"# Template: {template_name}",
        "# Review before submitting. This file is NOT auto-submitted by Ultimate.",
        "# Raw FASTQ inputs are read-only and are never copied or modified here.",
        *sbatch_lines,
        "",
        "set -euo pipefail",
        "",
        f"SAMPLESHEET={_quote(samplesheet)}",
        f"FASTQ_DIR={_quote(fastq_dir)}",
        f"INDEX={_quote(index)}",
        f"OUT_DIR={_quote(out_dir)}",
        f"THREADS=${{SLURM_CPUS_PER_TASK:-{threads}}}",
        "export THREADS",
        "",
    ]
    return "\n".join(preamble) + "\n" + "\n".join(rest)


def _quote(path: Path) -> str:
    text = str(path)
    if all(ch.isalnum() or ch in "/._-" for ch in text):
        return text
    escaped = text.replace('"', '\\"')
    return f'"{escaped}"'


def _render_expected_import_project(*, route: str, out_dir: Path, samplesheet: Path) -> str:
    if route == "salmon":
        body = [
            "project:",
            "  name: upstream_salmon_import",
            f"  output_dir: {out_dir}/runs/upstream_import",
            "  organism: human",
            "modules:",
            "  rnaseq:",
            "    enabled: true",
            f"    input_matrix: {out_dir}/counts_raw.tsv",
            "    analysis_level: smoke_backend",
        ]
    else:
        body = [
            "project:",
            "  name: upstream_starsolo_import",
            f"  output_dir: {out_dir}/runs/upstream_import",
            "  organism: human",
            "modules:",
            "  scrna:",
            "    enabled: true",
            f"    input_path: {out_dir}/star",
            "    input_type: 10x_mtx",
            "    analysis_level: smoke_backend",
        ]
    body.extend(
        [
            "# This is the *expected* import config after the upstream job succeeds.",
            f"# Source samplesheet (read-only): {samplesheet}",
            "upstream_package_note: generated_by_ultimate_prepare_upstream_review_before_submit",
        ]
    )
    return "\n".join(body) + "\n"


def prepare_upstream_package(
    *,
    route: str,
    samplesheet: Path,
    output_dir: Path,
    index: Path | None = None,
    fastq_dir: Path | None = None,
    threads: int | None = None,
    config: dict[str, Any] | None = None,
) -> PrepareUpstreamResult:
    """Build the job package; never submits anything."""
    route_key = str(route).strip().lower()
    output_dir = output_dir.expanduser().resolve()
    samplesheet = samplesheet.expanduser().resolve()
    fastq_base = (fastq_dir.expanduser().resolve() if fastq_dir else samplesheet.parent)
    index_path = index.expanduser().resolve() if index else None

    result = PrepareUpstreamResult(status="blocked", route=route_key, output_dir=output_dir)
    if route_key not in ROUTES:
        result.blocked_reason = f"unsupported_route:{route} (expected salmon or starsolo)"
        result.blockers.append(result.blocked_reason)
        return result

    budget = thread_budget(config or {}, threads)
    result.threads = budget.threads

    rows, samplesheet_error = _read_samplesheet(samplesheet)
    if samplesheet_error:
        result.blocked_reason = samplesheet_error
        result.blockers.append(samplesheet_error)
        return result
    result.n_samples = len(rows)

    if index_path is None:
        result.blocked_reason = f"missing_index:{ROUTES[route_key]['requires']}"
        result.blockers.append(result.blocked_reason)
        return result

    input_blockers = _validate_inputs(rows, fastq_base, index_path, route_key)
    raw_root = os.environ.get("ULTIMATE_RAW_ROOT", "")
    if raw_root:
        resolved_raw = Path(raw_root).resolve()
        if _is_inside(output_dir, resolved_raw):
            input_blockers.append("output_dir_inside_read_only_raw_tree")
    # FASTQ directories are inputs: writing the package inside the directory that
    # actually holds the reads would risk polluting a curated input tree even
    # when ULTIMATE_RAW_ROOT is unset. Only the resolved FASTQ parents count, so
    # an output directory beside the samplesheet stays allowed.
    fastq_parents = {
        _resolve_fastq(row[key], fastq_base).parent
        for row in rows
        for key in ("fastq_1", "fastq_2")
        if row.get(key) and row.get(key) != "-"
    }
    if any(_is_inside(output_dir, parent) for parent in fastq_parents):
        input_blockers.append("output_dir_inside_fastq_input_tree")
    if input_blockers:
        result.blocked_reason = input_blockers[0]
        result.blockers = input_blockers
        return result

    output_dir.mkdir(parents=True, exist_ok=True)
    template_name = ROUTES[route_key]["template"]
    try:
        template_path = resolve_packaged_asset("slurm", template_name)
    except FileNotFoundError:
        result.blocked_reason = f"missing_slurm_template:{template_name}"
        result.blockers.append(result.blocked_reason)
        return result
    template_body = template_path.read_text(encoding="utf-8")

    run_path = output_dir / "run_upstream.sbatch"
    project_path = output_dir / "expected_import_project.yaml"
    checksums_path = output_dir / "checksum_manifest.tsv"
    manifest_path = output_dir / "job_package_manifest.json"

    run_path.write_text(
        _render_run_upstream(
            route=route_key, samplesheet=samplesheet, fastq_dir=fastq_base, index=index_path,
            out_dir=output_dir / "upstream_out", threads=budget.threads,
            template_body=template_body, template_name=template_name,
        ),
        encoding="utf-8",
    )
    run_path.chmod(run_path.stat().st_mode | 0o111)
    project_path.write_text(
        _render_expected_import_project(route=route_key, out_dir=output_dir, samplesheet=samplesheet),
        encoding="utf-8",
    )

    checksum_rows = [("sample_id", "role", "path", "sha256", "read_only")]
    for row in rows:
        for key in ("fastq_1", "fastq_2"):
            value = row.get(key, "")
            if not value or value == "-":
                continue
            path = _resolve_fastq(value, fastq_base)
            checksum_rows.append((row["sample_id"], key, str(path), _checksum(path), "true"))
    with checksums_path.open("w", encoding="utf-8") as handle:
        for row in checksum_rows:
            handle.write("\t".join(str(item) for item in row) + "\n")

    result.status = "ready"
    result.artifacts = {
        "run_upstream_sbatch": str(run_path),
        "expected_import_project": str(project_path),
        "checksum_manifest": str(checksums_path),
        "job_package_manifest": str(manifest_path),
        "template": str(template_path),
    }
    manifest_path.write_text(
        json.dumps(
            {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "status": "ready",
                "route": route_key,
                "module": ROUTES[route_key]["module"],
                "samplesheet": str(samplesheet),
                "fastq_dir": str(fastq_base),
                "index": str(index_path),
                "output_dir": str(output_dir),
                "threads": budget.threads,
                "thread_source": budget.source,
                "raw_read_only": True,
                "auto_submitted": False,
                "submit_command": f"sbatch {run_path}",
                "artifacts": result.artifacts,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return result
