#!/usr/bin/env bash
#SBATCH --job-name=ult_scdna_backend
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" || ! -d "$PROJECT_ROOT/src/ultimate" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
export ULTIMATE_ROOT="$PROJECT_ROOT"
# Runtime artefacts live beside the code root on shared deployments; fall back to
# an in-checkout environment only when no runtime environment exists.
RUNTIME_ROOT="${ULTIMATE_RUNTIME_ROOT:-$(dirname -- "$PROJECT_ROOT")/ultimate_runtime}"
export ULTIMATE_RUNTIME_ROOT="$RUNTIME_ROOT"
if [[ -n "${ULTIMATE_CORE_ENV:-}" ]]; then
  CORE_ENV="$ULTIMATE_CORE_ENV"
else
  CORE_ENV=""
  for candidate in \
    "$RUNTIME_ROOT/envs/ultimate-core" \
    "$PROJECT_ROOT/.venv" \
    "$PROJECT_ROOT/.conda/envs/ultimate-core"
  do
    if [[ -x "$candidate/bin/python" ]]; then
      CORE_ENV="$candidate"
      break
    fi
  done
  if [[ -z "$CORE_ENV" ]]; then
    echo "No Ultimate core environment found. Build it with: sbatch slurm/xielab_setup_runtime_envs.sbatch" >&2
    exit 2
  fi
fi
PYTHON_BIN="${CORE_ENV}/bin/python"
OUT_DIR="${ULTIMATE_SCDNA_BACKEND_OUT:-$PROJECT_ROOT/validations/slurm_scdna_backend_0518}"
SOURCE_ROOT="${ULTIMATE_SCDNA_SOURCE:-${ULTIMATE_EXTERNAL_ROOT:-$PROJECT_ROOT/external}/0518}"

if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="python"
fi

write_skip_manifest() {
  mkdir -p "$OUT_DIR/reports"
  "$PYTHON_BIN" - "$OUT_DIR" "$1" "$2" <<'PY'
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

out = Path(sys.argv[1])
status = sys.argv[2]
reason = sys.argv[3]
manifest = {
    "generated_at": datetime.now(timezone.utc).isoformat(),
    "output_dir": str(out),
    "module": "scdna",
    "status": status,
    "analysis_level": "smoke_backend",
    "is_demo": False,
    "is_stub": True,
    "delivery_allowed": False,
    "validation_evidence_allowed": False,
    "non_delivery_reason": f"validation_not_completed:{reason}",
    "backend_id": "scdna.default.matrix_ready_handoff",
    "backend_status": "planned_fully_automatic",
    "backend_skip_reason": reason,
    "skip_reason": reason,
    "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
    "slurm_job_name": os.environ.get("SLURM_JOB_NAME", ""),
    "slurm": {
        "job_id": os.environ.get("SLURM_JOB_ID", ""),
        "job_name": os.environ.get("SLURM_JOB_NAME", ""),
        "submit_dir": os.environ.get("SLURM_SUBMIT_DIR", ""),
    },
}
(out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
(out / "reports" / "methods.md").write_text(
    "# scDNA backend validation skipped\n\n"
    "analysis_level: `smoke_backend`\n\n"
    "delivery_allowed: `false`\n\n"
    f"reason: {reason}\n",
    encoding="utf-8",
)
(out / "reports" / "report.html").write_text(
    f"<html><body><h1>scDNA backend validation skipped</h1><p>{reason}</p></body></html>",
    encoding="utf-8",
)
PY
}

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  write_skip_manifest "partial:dependency_required" "No Python runtime was found"
  exit 0
fi
if [[ ! -s "$SOURCE_ROOT/analysis_dna/stats/dna_mt_depth_summary.tsv" || ! -s "$SOURCE_ROOT/analysis_dna/stats/mapped_unmapped_summary.tsv" ]]; then
  write_skip_manifest "partial:data_required" "Missing 0518 scDNA/DNA validation tables under $SOURCE_ROOT"
  exit 0
fi

export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$PROJECT_ROOT/logs" "$OUT_DIR/config"

cat > "$OUT_DIR/config/project.yaml" <<YAML
project:
  name: slurm_scdna_backend_0518
  organism: human
  output_dir: $OUT_DIR
  server_root: $PROJECT_ROOT
  run_mode: interactive
  overwrite: true
samples:
  items:
    - sample_id: scdna_0518
      condition: internal_validation
      input_path: $SOURCE_ROOT
design:
  condition_column: condition
  control: internal_validation
  case: internal_validation
modules:
  scdna:
    enabled: true
    preset: standard
    source_root: $SOURCE_ROOT
    analysis_level: validated_backend
    validation_dataset: 0518 internal scDNA/RCA DNA QC results
    raw:
      enabled: false
      input_type: dna_qc_result_tables
report:
  title: scDNA unified backend validation
  style: soft_color
YAML

"$PYTHON_BIN" -m ultimate.cli preflight --config "$OUT_DIR/config/project.yaml" > "$OUT_DIR/preflight_stdout.json"
"$PYTHON_BIN" -m ultimate.cli run --config "$OUT_DIR/config/project.yaml" > "$OUT_DIR/run_stdout.json"
