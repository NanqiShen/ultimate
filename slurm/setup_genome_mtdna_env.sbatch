#!/usr/bin/env bash
#SBATCH --job-name=ult_env_gmtdna
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" || ! -d "$PROJECT_ROOT/src/ultimate" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
export ULTIMATE_ROOT="$PROJECT_ROOT"
bash "$PROJECT_ROOT/slurm/setup_singlecell_envs.sbatch" genome_mtdna

