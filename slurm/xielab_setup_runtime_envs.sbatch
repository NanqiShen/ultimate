#!/usr/bin/env bash
#SBATCH --job-name=ult_xielab_envs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=xielab
#SBATCH --partition=compute_pro
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00

# Build the three XieLab runtime environments with micromamba.
#
# This script only creates environments under $ULTIMATE_RUNTIME_ROOT/envs; it
# never submits jobs and never touches raw data. Login nodes are intentionally
# avoided because environment solving is a heavy operation.
#
# Usage:
#   ULTIMATE_ROOT=/dshare/xielab/analysisdata/TFC/shennq/ultimate \
#     sbatch slurm/xielab_setup_runtime_envs.sbatch [core] [scrna] [bio-r]
#   (no arguments builds all three)

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

echo "[xielab_envs] project_root=$PROJECT_ROOT"
echo "[xielab_envs] runtime_root=$RUNTIME_ROOT"
echo "[xielab_envs] env_root=$ENV_ROOT"

mkdir -p "$ENV_ROOT" "$RUNTIME_ROOT/logs" "$MICROMAMBA_ROOT/pkgs"

if [[ ! -x "$MAMBA" ]]; then
  echo "micromamba not found at $MAMBA" >&2
  echo "Install it under $RUNTIME_ROOT/micromamba first; nothing is downloaded automatically." >&2
  exit 2
fi

# Solving in parallel envs would fight over the package cache; keep it serial.
export MAMBA_ROOT_PREFIX="$MICROMAMBA_ROOT"
export CONDA_PKGS_DIRS="$MICROMAMBA_ROOT/pkgs"
export MAMBA_DOWNLOAD_THREADS="${SLURM_CPUS_PER_TASK:-2}"

declare -A ENV_FILE=(
  [core]="envs/xielab.core.yml"
  [scrna]="envs/xielab.scrna.py.yml"
  [bio-r]="envs/xielab.bio.r.yml"
)
declare -A ENV_NAME=(
  [core]="ultimate-core"
  [scrna]="ultimate-scrna"
  [bio-r]="ultimate-bio-r"
)

requested=("$@")
if [[ ${#requested[@]} -eq 0 ]]; then
  requested=(core scrna bio-r)
fi

status_file="$RUNTIME_ROOT/manifests/xielab_env_build_${SLURM_JOB_ID:-manual}.tsv"
printf 'env\tstatus\tprefix\tlog\n' > "$status_file"

overall=0
for key in "${requested[@]}"; do
  file="${ENV_FILE[$key]:-}"
  name="${ENV_NAME[$key]:-}"
  if [[ -z "$file" || -z "$name" ]]; then
    echo "[xielab_envs] unknown environment selector: $key (expected core, scrna, bio-r)" >&2
    printf '%s\tskipped\t\tunknown_selector\n' "$key" >> "$status_file"
    overall=2
    continue
  fi
  definition="$PROJECT_ROOT/$file"
  if [[ ! -f "$definition" ]]; then
    echo "[xielab_envs] definition missing: $definition" >&2
    printf '%s\tskipped\t\tmissing_definition:%s\n' "$key" "$definition" >> "$status_file"
    overall=2
    continue
  fi
  prefix="$ENV_ROOT/$name"
  log="$RUNTIME_ROOT/logs/xielab_env_${name}_${SLURM_JOB_ID:-manual}.log"
  echo "[xielab_envs] building $name from $file -> $prefix (log: $log)"
  if "$MAMBA" env update --prefix "$prefix" --file "$definition" --prune > "$log" 2>&1 \
     || "$MAMBA" create --prefix "$prefix" --file "$definition" --yes >> "$log" 2>&1; then
    printf '%s\tready\t%s\t%s\n' "$key" "$prefix" "$log" >> "$status_file"
  else
    printf '%s\tfailed\t%s\t%s\n' "$key" "$prefix" "$log" >> "$status_file"
    echo "[xielab_envs] FAILED building $name; see $log" >&2
    overall=1
  fi
done

echo "[xielab_envs] status table: $status_file"
cat "$status_file"
echo "[xielab_envs] verify with:"
echo "  $MAMBA run -p $ENV_ROOT/ultimate-core python -c 'import ultimate'"
echo "  $ENV_ROOT/ultimate-bio-r/bin/Rscript -e 'library(limma)'"
exit "$overall"
