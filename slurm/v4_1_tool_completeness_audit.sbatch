#!/usr/bin/env bash
#SBATCH --job-name=ult_v4_1_tools
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00

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

mkdir -p "$PROJECT_ROOT/logs" "$PROJECT_ROOT/audits" "$PROJECT_ROOT/reports"
cd "$PROJECT_ROOT"
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

"$CORE_ENV/bin/python" -m ultimate.cli audit-tools \
  --root "$PROJECT_ROOT" \
  --output-dir "$PROJECT_ROOT/audits/tools_v4_1_latest"

"$CORE_ENV/bin/python" -m ultimate.cli validation-index \
  --root "$PROJECT_ROOT" \
  --output-dir "$PROJECT_ROOT/reports/validation_index_v4_1"

"$CORE_ENV/bin/python" -m ultimate.cli audit-production \
  --root "$PROJECT_ROOT" \
  --output-dir "$PROJECT_ROOT/audits/production_v4_1_latest"

"$CORE_ENV/bin/python" -m ultimate.cli audit-backends \
  --root "$PROJECT_ROOT" \
  --output-dir "$PROJECT_ROOT/audits/backends_v4_1_latest"

"$CORE_ENV/bin/python" -m ultimate.cli tool-completeness \
  --root "$PROJECT_ROOT" \
  --output-dir "$PROJECT_ROOT/audits/tool_completeness_latest"

"$CORE_ENV/bin/python" -m ultimate.cli order-readiness \
  --root "$PROJECT_ROOT" \
  --output-dir "$PROJECT_ROOT/audits/order_readiness_latest"

if [[ -f "$PROJECT_ROOT/01_tools/storage_audit.py" ]]; then
  "$CORE_ENV/bin/python" "$PROJECT_ROOT/01_tools/storage_audit.py" \
    --root "$PROJECT_ROOT" \
    --output-dir "$PROJECT_ROOT/audits/storage_v4_1_latest" \
    --budget-gb 500
fi

cp "$PROJECT_ROOT/audits/tool_completeness_latest/v4_1_tool_completeness_report.md" \
  "$PROJECT_ROOT/reports/v4_1_tool_completeness_report.md"

