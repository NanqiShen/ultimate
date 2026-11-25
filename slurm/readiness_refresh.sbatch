#!/usr/bin/env bash
#SBATCH --job-name=ult_ready_refresh
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

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
VALIDATION_INDEX_OUT="${ULTIMATE_VALIDATION_INDEX_OUTPUT:-$PROJECT_ROOT/reports/validation_index}"
AUDIT_OUT="${ULTIMATE_PRODUCTION_AUDIT_OUTPUT:-$PROJECT_ROOT/audits/production_readiness_latest}"
GUARD_OUT="${ULTIMATE_VALIDATION_GUARD_OUTPUT:-$PROJECT_ROOT/audits/validation_guard_latest/validation_manifest_guard.tsv}"
METHODS_OUT="${ULTIMATE_VALIDATION_METHODS_OUTPUT:-$PROJECT_ROOT/audits/validation_guard_latest/validation_methods_backfill.tsv}"
MVP_OUT="${ULTIMATE_VALIDATION_MVP_OUTPUT:-$PROJECT_ROOT/audits/validation_guard_latest/validation_mvp_artifacts.tsv}"

mkdir -p "$PROJECT_ROOT/logs" "$VALIDATION_INDEX_OUT" "$AUDIT_OUT" "$(dirname "$GUARD_OUT")"

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "Missing ultimate-core Python: $PYTHON_BIN" >&2
  exit 2
fi

export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

"$PYTHON_BIN" "$PROJECT_ROOT/01_tools/ensure_validation_methods.py" \
  --root "$PROJECT_ROOT" \
  --validations-dir "$PROJECT_ROOT/validations" \
  --output-tsv "$METHODS_OUT"

"$PYTHON_BIN" "$PROJECT_ROOT/01_tools/ensure_validation_mvp_artifacts.py" \
  --root "$PROJECT_ROOT" \
  --validations-dir "$PROJECT_ROOT/validations" \
  --output-tsv "$MVP_OUT"

"$PYTHON_BIN" -m ultimate.cli validation-index \
  --root "$PROJECT_ROOT" \
  --output-dir "$VALIDATION_INDEX_OUT"

"$PYTHON_BIN" -m ultimate.cli audit-production \
  --root "$PROJECT_ROOT" \
  --output-dir "$AUDIT_OUT"

"$PYTHON_BIN" "$PROJECT_ROOT/01_tools/check_validation_manifests.py" \
  --root "$PROJECT_ROOT" \
  --validations-dir "$PROJECT_ROOT/validations" \
  --output-tsv "$GUARD_OUT"
