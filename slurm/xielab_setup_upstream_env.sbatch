#!/usr/bin/env bash
#SBATCH --job-name=ult_xielab_upstream
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=xielab
#SBATCH --partition=compute_pro
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00

# Build the XieLab open-source upstream quantification environment.
#
# Creates $ULTIMATE_RUNTIME_ROOT/envs/ultimate-upstream from envs/xielab.upstream.yml
# (fastqc/multiqc/fastp/salmon/star/subread/samtools/alevin-fry). This replaces the
# need for nextflow/apptainer/cellranger, which are not installed on compute_pro.
#
# This script only creates an environment: it never submits jobs, never writes to
# /dshare/dhome and never touches raw FASTQ under /dshare/xielab/rawdata.
#
# Usage:
#   ULTIMATE_ROOT=/dshare/xielab/analysisdata/TFC/shennq/ultimate \
#     sbatch slurm/xielab_setup_upstream_env.sbatch

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
RUNTIME_ROOT="${ULTIMATE_RUNTIME_ROOT:-$(dirname -- "$PROJECT_ROOT")/ultimate_runtime}"
ENV_ROOT="${ULTIMATE_ENV_ROOT:-$RUNTIME_ROOT/envs}"
MICROMAMBA_ROOT="${ULTIMATE_MICROMAMBA_ROOT:-$RUNTIME_ROOT/micromamba-root}"
MAMBA="$RUNTIME_ROOT/micromamba/bin/micromamba"
ENV_PREFIX="$ENV_ROOT/ultimate-upstream"
DEFINITION="$PROJECT_ROOT/envs/xielab.upstream.yml"
LOG_DIR="$RUNTIME_ROOT/logs"
LOG="$LOG_DIR/xielab_upstream_env_${SLURM_JOB_ID:-manual}.log"
STATUS_TABLE="$RUNTIME_ROOT/manifests/xielab_upstream_env_${SLURM_JOB_ID:-manual}.tsv"

echo "[xielab_upstream] project_root=$PROJECT_ROOT"
echo "[xielab_upstream] env_prefix=$ENV_PREFIX"

mkdir -p "$ENV_ROOT" "$LOG_DIR" "$MICROMAMBA_ROOT/pkgs" "$RUNTIME_ROOT/manifests"

if [[ ! -f "$DEFINITION" ]]; then
  echo "[xielab_upstream] definition missing: $DEFINITION" >&2
  exit 2
fi
if [[ ! -x "$MAMBA" ]]; then
  echo "[xielab_upstream] micromamba not found at $MAMBA" >&2
  echo "Install micromamba under $RUNTIME_ROOT/micromamba first; nothing is downloaded automatically." >&2
  exit 2
fi

export MAMBA_ROOT_PREFIX="$MICROMAMBA_ROOT"
export CONDA_PKGS_DIRS="$MICROMAMBA_ROOT/pkgs"
export MAMBA_DOWNLOAD_THREADS="${SLURM_CPUS_PER_TASK:-2}"

printf 'env\tstatus\tprefix\tlog\n' > "$STATUS_TABLE"

echo "[xielab_upstream] building ultimate-upstream from envs/xielab.upstream.yml (log: $LOG)"
status="failed"
if "$MAMBA" env update --prefix "$ENV_PREFIX" --file "$DEFINITION" --prune > "$LOG" 2>&1 \
   || "$MAMBA" create --prefix "$ENV_PREFIX" --file "$DEFINITION" --yes >> "$LOG" 2>&1; then
  status="ready"
fi
printf 'ultimate-upstream\t%s\t%s\t%s\n' "$status" "$ENV_PREFIX" "$LOG" >> "$STATUS_TABLE"
cat "$STATUS_TABLE"

if [[ "$status" != "ready" ]]; then
  echo "[xielab_upstream] FAILED; see $LOG" >&2
  exit 1
fi

echo "[xielab_upstream] verifying tool availability (threads bounded by SLURM_CPUS_PER_TASK)"
for tool in fastqc multiqc fastp salmon STAR featureCounts samtools alevin-fry; do
  if "$ENV_PREFIX/bin/$tool" --version > /dev/null 2>&1 || command -v "$tool" > /dev/null 2>&1; then
    echo "  OK   $tool"
  else
    echo "  WARN $tool not resolved on PATH inside the environment"
  fi
done
echo "[xielab_upstream] done. Activate with:"
echo "  export ULTIMATE_ENV=$ENV_PREFIX"
