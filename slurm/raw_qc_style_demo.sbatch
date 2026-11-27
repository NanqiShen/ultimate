#!/usr/bin/env bash
#SBATCH --job-name=ult_raw_qc_style
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" || ! -d "$PROJECT_ROOT/src/ultimate" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
export ULTIMATE_ROOT="$PROJECT_ROOT"
CONFIG_PATH="${ULTIMATE_CONFIG:-$PROJECT_ROOT/projects/raw_qc_style_demo/config/project.yaml}"
ENV_PREFIX="${ULTIMATE_ENV:-$PROJECT_ROOT/.conda/envs/ultimate-core}"

mkdir -p "$PROJECT_ROOT/logs"
source "${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
conda activate "$ENV_PREFIX"
cd "$PROJECT_ROOT"
ultimate run --config "$CONFIG_PATH"
