#!/usr/bin/env bash
#SBATCH --job-name=ult_val_index
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
CORE_ENV="${ULTIMATE_CORE_ENV:-$PROJECT_ROOT/.conda/envs/ultimate-core}"
OUTPUT_DIR="${ULTIMATE_VALIDATION_INDEX_OUTPUT:-$PROJECT_ROOT/reports/validation_index}"

mkdir -p "$PROJECT_ROOT/logs" "$OUTPUT_DIR"
cd "$PROJECT_ROOT"

PYTHONPATH="$PROJECT_ROOT/src" "$CORE_ENV/bin/python" -m ultimate.cli validation-index --root "$PROJECT_ROOT" --output-dir "$OUTPUT_DIR"
