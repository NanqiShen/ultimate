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

# Sub-directory of the package that the Slurm template writes artefacts into.
ARTIFACTS_SUBDIR = "upstream_out"

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


# Accepted column aliases (case-insensitive, matched after trimming).
SAMPLE_ALIASES = ("sample", "sample_id", "sampleid", "sample name", "sample_name")
FASTQ1_ALIASES = ("fastq_1", "fq1", "r1", "fastq_r1", "read1", "fastq1")
FASTQ2_ALIASES = ("fastq_2", "fq2", "r2", "fastq_r2", "read2", "fastq2")
STRANDEDNESS_ALIASES = ("strandedness", "strand", "orientation")
CONDITION_ALIASES = ("condition", "group", "treatment", "phenotype", "status")

NORMALIZED_COLUMNS = ("sample_id", "fastq_1", "fastq_2", "strandedness", "condition")


def _pick(clean: dict[str, str], aliases: tuple[str, ...]) -> str:
    """First non-empty value among ``aliases`` (column names are case-insensitive)."""
    for alias in aliases:
        value = clean.get(alias, "")
        if value:
            return value
    return ""


def _normalize_columns(row: dict[str, str]) -> dict[str, str]:
    clean: dict[str, str] = {}
    for key, value in row.items():
        if key is None:
            continue
        clean[str(key).strip().lower().replace("-", "_")] = str(value or "").strip()
    return clean


def _read_samplesheet(path: Path) -> tuple[list[dict[str, str]], str]:
    """Parse and normalise the FASTQ samplesheet; return (rows, error).

    The same normalised rows are written to ``samplesheet_normalized.tsv`` and are
    the *only* thing the generated Slurm job consumes, so the validator and the
    job script can never disagree about column order or aliases.
    """
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
    seen: set[str] = set()
    for index, row in enumerate(rows, start=2):
        clean = _normalize_columns(row)
        # Skip fully blank lines rather than failing on them.
        if not any(clean.values()):
            continue
        sample = _pick(clean, SAMPLE_ALIASES)
        fq1 = _pick(clean, FASTQ1_ALIASES)
        if not sample:
            return [], f"samplesheet_missing_sample_at_row:{index}"
        if not fq1:
            return [], f"samplesheet_missing_fastq_1_at_row:{index}"
        if sample in seen:
            return [], f"duplicate_sample_id:{sample}"
        seen.add(sample)
        normalized.append(
            {
                "sample_id": sample,
                "fastq_1": fq1,
                "fastq_2": _pick(clean, FASTQ2_ALIASES),
                "strandedness": _pick(clean, STRANDEDNESS_ALIASES) or "auto",
                "condition": _pick(clean, CONDITION_ALIASES) or "control",
            }
        )
    if not normalized:
        return [], "empty_samplesheet:no_data_rows"
    return normalized, ""


def _validate_pairs(rows: list[dict[str, str]]) -> list[str]:
    """A sample either has both reads or neither; half-pairs are ambiguous."""
    blockers: list[str] = []
    for row in rows:
        has1 = bool(row.get("fastq_1") and row["fastq_1"] != "-")
        has2 = bool(row.get("fastq_2") and row["fastq_2"] != "-")
        if has2 and not has1:
            blockers.append(f"missing_paired_read:fastq_1:{row['sample_id']}")
        if has1 and has2:
            continue
        # Single-end is legitimate; nothing to add.
    return blockers


def _write_normalized_samplesheet(path: Path, rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(NORMALIZED_COLUMNS)
        for row in rows:
            writer.writerow([row.get(col, "") for col in NORMALIZED_COLUMNS])
    return path


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


def _count_semantics(route: str, artifacts_dir: Path) -> dict[str, str]:
    """State exactly what the produced counts are - and what they are not."""
    if route == "salmon":
        return {
            "level": "transcript",
            "unit": "salmon NumReads rounded to integers",
            "is_gene_level_raw_integer_counts": "false",
            "raw_estimates_kept_at": str(artifacts_dir / "quant"),
            "note": (
                "NumReads is a transcript-level abundance estimate, not gene-level raw "
                "integer counts. Preserve quant.sf and map transcripts to genes before any "
                "gene-level differential expression."
            ),
        }
    return {
        "level": "cell_x_gene",
        "unit": "STARsolo UMI counts in the standard 10x filtered/ matrix",
        "is_gene_level_raw_integer_counts": "true",
        "raw_estimates_kept_at": str(artifacts_dir / "star"),
        "note": (
            "STAR --soloType CB_UMI_Simple emits non-negative integer UMI counts; the "
            "matrix is the standard 10x matrix.mtx format."
        ),
    }


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


def _render_expected_import_project(
    *,
    route: str,
    artifacts_dir: Path,
    samplesheet: Path,
    sample_ids: list[str],
) -> str:
    """Render the import config that consumes the job's *actual* output paths.

    ``artifacts_dir`` is the directory the Slurm template writes into
    (``<package>/upstream_out``); using the package root here would point the
    import at paths that never exist.
    """
    run_root = artifacts_dir / "runs" / "upstream_import"
    if route == "salmon":
        body = [
            "project:",
            "  name: upstream_salmon_import",
            f"  output_dir: {run_root}",
            "  organism: human",
            "samples:",
            f"  samplesheet: {samplesheet}",
            "modules:",
            "  rnaseq:",
            "    enabled: true",
            f"    input_matrix: {artifacts_dir}/counts_raw.tsv",
            f"    samplesheet: {samplesheet}",
            "    analysis_level: smoke_backend",
            "# Count semantics (do not overstate):",
            "#   counts_raw.tsv holds salmon transcript-level NumReads values rounded to",
            "#   integers. They are NOT gene-level raw integer counts; map transcripts to",
            "#   genes before any gene-level differential expression. The untouched",
            "#   per-sample quant.sf files are kept under quant/<sample>/quant.sf.",
            f"# upstream_raw_quant: {artifacts_dir}/quant",
            f"# upstream_import_config: {artifacts_dir}/import_config.yaml",
        ]
    else:
        # STARsolo writes one filtered/ 10x directory per sample; the import must
        # point at that directory, not at its parent.
        per_sample = "\n".join(
            f"    - sample_id: {sample}\n      input_path: {artifacts_dir}/star/{sample}/filtered"
            for sample in sample_ids
        )
        body = [
            "project:",
            "  name: upstream_starsolo_import",
            f"  output_dir: {run_root}",
            "  organism: human",
            "modules:",
            "  scrna:",
            "    enabled: true",
            f"    input_path: {artifacts_dir}/star/{sample_ids[0]}/filtered",
            "    input_type: 10x_mtx",
            "    analysis_level: smoke_backend",
            "    samples:",
            per_sample,
            "# Each input_path above is the standard 10x filtered/ directory produced",
            "# by STARsolo for that sample (matrix.mtx.gz / barcodes.tsv.gz /",
            "# features.tsv.gz).",
            f"# upstream_import_config: {artifacts_dir}/import_config.yaml",
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
    input_blockers.extend(_validate_pairs(rows))
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

    # Never silently overwrite a previous run: an existing non-empty package is
    # treated as a historical result and must be moved or removed by a human.
    if output_dir.is_dir() and any(output_dir.iterdir()):
        existing = sorted(p.name for p in output_dir.iterdir())
        result.blocked_reason = f"output_dir_not_empty:{','.join(existing[:5])}"
        result.blockers = [result.blocked_reason]
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

    artifacts_dir = output_dir / ARTIFACTS_SUBDIR
    # Single source of truth: the job reads only this normalised table, so the
    # validator and the sbatch template can never disagree about columns.
    normalized_path = _write_normalized_samplesheet(output_dir / "samplesheet_normalized.tsv", rows)
    run_path.write_text(
        _render_run_upstream(
            route=route_key, samplesheet=normalized_path, fastq_dir=fastq_base, index=index_path,
            out_dir=artifacts_dir, threads=budget.threads,
            template_body=template_body, template_name=template_name,
        ),
        encoding="utf-8",
    )
    run_path.chmod(run_path.stat().st_mode | 0o111)
    project_path.write_text(
        _render_expected_import_project(
            route=route_key,
            artifacts_dir=artifacts_dir,
            samplesheet=normalized_path,
            sample_ids=[row["sample_id"] for row in rows],
        ),
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
        "samplesheet_normalized": str(normalized_path),
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
                "samplesheet_source": str(samplesheet),
                "samplesheet_normalized": str(normalized_path),
                "fastq_dir": str(fastq_base),
                "index": str(index_path),
                "output_dir": str(output_dir),
                "threads": budget.threads,
                "thread_source": budget.source,
                "raw_read_only": True,
                "auto_submitted": False,
                "count_semantics": _count_semantics(route_key, artifacts_dir),
                "submit_command": f"sbatch {run_path}",
                "artifacts": result.artifacts,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return result
