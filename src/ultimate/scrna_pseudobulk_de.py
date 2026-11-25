from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any


BACKEND_ID = "scrna.pseudobulk.edger_paired"


def run_scrna_pseudobulk_de(
    *,
    output_dir: Path,
    sce_rda: Path | None = None,
    counts: Path | None = None,
    design: Path | None = None,
    rscript: Path | None = None,
    sce_object: str | None = None,
    assay: str = "counts",
    donor_column: str = "ind",
    condition_column: str = "stim",
    cell_type_column: str = "cell",
    cell_type: str,
    control: str,
    case: str,
    multiplet_column: str | None = None,
    multiplet_keep: str | None = None,
    min_cells: int = 1,
    min_pairs: int = 4,
    analysis_level: str = "validated_backend",
    robust: bool = True,
) -> dict[str, Any]:
    """Run the paired edgeR pseudobulk DE validation backend through Ultimate."""

    if sce_rda is None and (counts is None or design is None):
        raise ValueError("Provide either --sce-rda or both --counts and --design.")
    if sce_rda is not None and (counts is not None or design is not None):
        raise ValueError("--sce-rda cannot be combined with --counts/--design.")
    if min_cells < 1:
        raise ValueError("--min-cells must be at least 1.")
    if min_pairs < 2:
        raise ValueError("--min-pairs must be at least 2.")

    output_dir = output_dir.resolve()
    tables_dir = output_dir / "results" / "tables"
    figures_dir = output_dir / "results" / "figures"
    objects_dir = output_dir / "objects"
    logs_dir = output_dir / "logs"
    for directory in (tables_dir, figures_dir, objects_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    rscript_found = shutil.which("Rscript")
    if rscript is None and not rscript_found:
        raise RuntimeError("Rscript was not found on PATH; provide --rscript.")
    rscript_path = rscript or Path(rscript_found)
    backend_script = Path(__file__).resolve().parents[2] / "scripts" / "R" / "scrna_paired_pseudobulk_edger.R"
    if not backend_script.exists():
        raise RuntimeError(f"Backend script is missing: {backend_script}")

    command = [
        str(rscript_path),
        str(backend_script),
        "--tables-dir",
        str(tables_dir),
        "--backend-id",
        BACKEND_ID,
        "--analysis-level",
        analysis_level,
        "--donor-column",
        donor_column,
        "--condition-column",
        condition_column,
        "--cell-type-column",
        cell_type_column,
        "--cell-type",
        cell_type,
        "--control",
        control,
        "--case",
        case,
        "--min-cells",
        str(min_cells),
        "--min-pairs",
        str(min_pairs),
        "--robust",
        "true" if robust else "false",
    ]
    if sce_rda is not None:
        command.extend(["--sce-rda", str(sce_rda.resolve()), "--assay", assay])
        if sce_object:
            command.extend(["--sce-object", sce_object])
        if multiplet_column and multiplet_keep:
            command.extend(["--multiplet-column", multiplet_column, "--multiplet-keep", multiplet_keep])
    else:
        command.extend(["--counts", str(counts.resolve()), "--design", str(design.resolve())])

    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    log_path = logs_dir / "scrna_pseudobulk_de_backend.log"
    log_path.write_text(
        "COMMAND\t" + " ".join(command) + "\n\nSTDOUT\n" + completed.stdout + "\n\nSTDERR\n" + completed.stderr,
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError(f"scrna_pseudobulk_de_backend_failed:{completed.returncode}:{log_path}")
    manifest_path = tables_dir / "pseudobulk_de_backend_manifest.json"
    if not manifest_path.exists():
        raise RuntimeError(f"pseudobulk_de_backend_manifest_missing:{manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["ultimate_cli_command"] = command
    manifest["ultimate_cli_log"] = str(log_path)
    manifest["output_dir"] = str(output_dir)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest
