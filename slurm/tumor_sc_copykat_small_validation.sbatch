#!/usr/bin/env bash
#SBATCH --job-name=ult_copykat_small
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" || ! -d "$PROJECT_ROOT/src/ultimate" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
export ULTIMATE_ROOT="$PROJECT_ROOT"
export ULTIMATE_TUMOR_SC_INPUT="${ULTIMATE_TUMOR_SC_INPUT:-${ULTIMATE_EXTERNAL_ROOT:-$PROJECT_ROOT/external}/nsclc_virtual_sc/02_processed_scRNA/Maynard_raw.h5ad}"
export ULTIMATE_TUMOR_SC_OUT="${ULTIMATE_TUMOR_SC_OUT:-$PROJECT_ROOT/validations/slurm_tumor_sc_copykat_small}"
export ULTIMATE_TUMOR_SC_MAX_CELLS="${ULTIMATE_TUMOR_SC_MAX_CELLS:-3000}"
export ULTIMATE_TUMOR_SC_RUN_COPYKAT="${ULTIMATE_TUMOR_SC_RUN_COPYKAT:-1}"
export ULTIMATE_TUMOR_SC_RUN_INFERCNV="${ULTIMATE_TUMOR_SC_RUN_INFERCNV:-1}"
export ULTIMATE_TUMOR_SC_APPROVE_FULL_CNV="${ULTIMATE_TUMOR_SC_APPROVE_FULL_CNV:-1}"
export ULTIMATE_TUMOR_SC_COPYKAT_MAX_CELLS="${ULTIMATE_TUMOR_SC_COPYKAT_MAX_CELLS:-500}"
export ULTIMATE_TUMOR_SC_COPYKAT_MAX_GENES="${ULTIMATE_TUMOR_SC_COPYKAT_MAX_GENES:-6000}"
export ULTIMATE_TUMOR_SC_GENE_ORDER_REF="${ULTIMATE_TUMOR_SC_GENE_ORDER_REF:-${ULTIMATE_EXTERNAL_ROOT:-$PROJECT_ROOT/external}/nsclc_virtual_sc/02_processed_scRNA/Lambrechts_all_harmonized.h5ad}"

bash "$PROJECT_ROOT/slurm/tumor_sc_validation.sbatch"
