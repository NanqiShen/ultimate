#!/usr/bin/env bash
#SBATCH --job-name=ult_gapfill_validate
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
PYTHON="${CORE_ENV}/bin/python"
if [[ ! -x "$PYTHON" ]]; then
  PYTHON="python"
fi

mkdir -p "$PROJECT_ROOT/logs" "$PROJECT_ROOT/validations" "$PROJECT_ROOT/audits"
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

"$PYTHON" "$PROJECT_ROOT/01_tools/validate_perturb_seq_demo.py" \
  --output-dir "$PROJECT_ROOT/validations/slurm_perturb_seq_demo"

PERTURB_PUBLIC_H5AD="${ULTIMATE_PERTURB_SEQ_PUBLIC_H5AD:-$PROJECT_ROOT/public_data/perturb_seq/adamson_2016_pilot/adamson_2016_pilot.h5ad}"
if [[ -s "$PERTURB_PUBLIC_H5AD" ]]; then
  PERTURB_PYTHON="${ULTIMATE_SCRNA_ENV:-$PROJECT_ROOT/.conda/envs/ultimate-scrna}/bin/python"
  if [[ ! -x "$PERTURB_PYTHON" ]]; then
    PERTURB_PYTHON="$PYTHON"
  fi
  "$PERTURB_PYTHON" "$PROJECT_ROOT/01_tools/validate_perturb_seq_demo.py" \
    --input-h5ad "$PERTURB_PUBLIC_H5AD" \
    --source-url "https://exampledata.scverse.org/pertpy/adamson_2016_pilot.h5ad" \
    --output-dir "$PROJECT_ROOT/validations/slurm_perturb_seq_adamson_public"
fi

"$PYTHON" "$PROJECT_ROOT/01_tools/validate_hto_demux_demo.py" \
  --output-dir "$PROJECT_ROOT/validations/slurm_hto_demux_demo"

HTO_PUBLIC_TABLE="${ULTIMATE_HTO_DEMUX_PUBLIC_TABLE:-$PROJECT_ROOT/public_data/hto_demux/seurat_hto12/hto12_hto_mtx_5000.tsv}"
if [[ -s "$HTO_PUBLIC_TABLE" ]]; then
  "$PYTHON" "$PROJECT_ROOT/01_tools/validate_hto_demux_demo.py" \
    --input-table "$HTO_PUBLIC_TABLE" \
    --source-url "https://satijalab.org/seurat/articles/hashing_vignette.html" \
    --output-dir "$PROJECT_ROOT/validations/slurm_hto_demux_seurat_public"
fi

"$PYTHON" "$PROJECT_ROOT/01_tools/validate_genotype_demux_demo.py" \
  --output-dir "$PROJECT_ROOT/validations/slurm_genotype_demux_demo"

GENOTYPE_PUBLIC_INPUT="${ULTIMATE_GENOTYPE_DEMUX_PUBLIC_INPUT:-$PROJECT_ROOT/public_data/genotype_demux/vireo/data/cellSNP_mat}"
if [[ -s "$GENOTYPE_PUBLIC_INPUT/cellSNP.tag.AD.mtx" && -s "$GENOTYPE_PUBLIC_INPUT/cellSNP.tag.DP.mtx" ]]; then
  "$PYTHON" "$PROJECT_ROOT/01_tools/validate_genotype_demux_demo.py" \
    --input-dir "$GENOTYPE_PUBLIC_INPUT" \
    --source-url "https://github.com/single-cell-genetics/vireo/tree/master/data/cellSNP_mat" \
    --output-dir "$PROJECT_ROOT/validations/slurm_genotype_demux_vireo_public"
fi

"$PYTHON" -m ultimate.cli audit-production \
  --root "$PROJECT_ROOT" \
  --output-dir "$PROJECT_ROOT/audits/gapfill_specialty_post_validation"
