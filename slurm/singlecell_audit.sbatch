#!/usr/bin/env bash
#SBATCH --job-name=ult_sc_audit
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" || ! -d "$PROJECT_ROOT/src/ultimate" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
export ULTIMATE_ROOT="$PROJECT_ROOT"
CORE_ENV="${ULTIMATE_CORE_ENV:-$PROJECT_ROOT/.conda/envs/ultimate-core}"
CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
OUTPUT_DIR="${ULTIMATE_AUDIT_OUTPUT:-$PROJECT_ROOT/audits/singlecell_slurm_latest}"

mkdir -p "$PROJECT_ROOT/logs" "$OUTPUT_DIR"
source "$CONDA_SH"
conda activate "$CORE_ENV"
cd "$PROJECT_ROOT"

ultimate audit-singlecell --root "$PROJECT_ROOT" --output-dir "$OUTPUT_DIR"

