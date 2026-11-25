#!/usr/bin/env bash
#SBATCH --job-name=ultimate
#SBATCH --output=ultimate_%j.out
#SBATCH --error=ultimate_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00

set -euo pipefail

CONFIG_PATH="${1:?Usage: ultimate_run.sbatch <config.yaml> [production_approval.json]}"
PRODUCTION_APPROVAL="${2:-}"
PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" || ! -d "$PROJECT_ROOT/src/ultimate" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
export ULTIMATE_ROOT="$PROJECT_ROOT"
CONFIG_PATH="$(readlink -f "$CONFIG_PATH")"

# Runtime artefacts live outside the code root on shared deployments
# (<site>/ultimate + <site>/ultimate_runtime).
RUNTIME_ROOT="${ULTIMATE_RUNTIME_ROOT:-$(dirname -- "$PROJECT_ROOT")/ultimate_runtime}"
export ULTIMATE_RUNTIME_ROOT="$RUNTIME_ROOT"
export ULTIMATE_ENV_ROOT="${ULTIMATE_ENV_ROOT:-$RUNTIME_ROOT/envs}"
export ULTIMATE_CACHE_ROOT="${ULTIMATE_CACHE_ROOT:-$RUNTIME_ROOT/cache}"

# First existing environment wins; ULTIMATE_ENV always overrides discovery.
discover_env_prefix() {
  if [[ -n "${ULTIMATE_ENV:-}" ]]; then
    [[ -x "$ULTIMATE_ENV/bin/python" ]] && { echo "$ULTIMATE_ENV"; return 0; }
    echo "ULTIMATE_ENV does not contain bin/python: $ULTIMATE_ENV" >&2
    return 1
  fi
  local candidate
  for candidate in \
    "$ULTIMATE_ENV_ROOT/ultimate-core" \
    "$PROJECT_ROOT/.venv" \
    "$PROJECT_ROOT/.conda/envs/ultimate-core"
  do
    if [[ -x "$candidate/bin/python" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

ENV_PREFIX="$(discover_env_prefix)" || {
  echo "No usable Ultimate environment found. Searched:" >&2
  echo "  ULTIMATE_ENV, $ULTIMATE_ENV_ROOT/ultimate-core, $PROJECT_ROOT/.venv, $PROJECT_ROOT/.conda/envs/ultimate-core" >&2
  echo "Build one with: sbatch slurm/xielab_setup_runtime_envs.sbatch" >&2
  exit 2
}
ENV_KIND="path_prefix"
CONDA_SH="${CONDA_SH:-}"

mkdir -p "$PROJECT_ROOT/logs"
CONFIG_DIR="$(dirname "$CONFIG_PATH")"
if [[ "$(basename "$CONFIG_DIR")" == "config" ]]; then
  JOB_DIR="$(dirname "$CONFIG_DIR")"
else
  JOB_DIR=""
fi
if [[ -n "$JOB_DIR" && "$JOB_DIR" == "$PROJECT_ROOT"/jobs/* ]]; then
  JOB_LOG_DIR="$JOB_DIR/logs"
else
  JOB_LOG_DIR="$PROJECT_ROOT/logs"
fi
mkdir -p "$JOB_LOG_DIR"
JOB_LOG_PREFIX="${SLURM_JOB_ID:-manual}_ultimate"
exec > >(tee -a "$JOB_LOG_DIR/${JOB_LOG_PREFIX}.out") 2> >(tee -a "$JOB_LOG_DIR/${JOB_LOG_PREFIX}.err" >&2)
cat > "$JOB_LOG_DIR/${JOB_LOG_PREFIX}_slurm_context.json" <<JSON
{
  "slurm_job_id": "${SLURM_JOB_ID:-}",
  "slurm_job_name": "${SLURM_JOB_NAME:-}",
  "slurm_submit_dir": "${SLURM_SUBMIT_DIR:-}",
  "slurm_cpus_per_task": "${SLURM_CPUS_PER_TASK:-}",
  "slurm_job_nodelist": "${SLURM_JOB_NODELIST:-}",
  "config_path": "$CONFIG_PATH",
  "job_dir": "$JOB_DIR",
  "project_root": "$PROJECT_ROOT"
}
JSON
MICROMAMBA_BIN="$RUNTIME_ROOT/micromamba/bin/micromamba"

# Activation must not depend on conda being on PATH: Xielab only ships micromamba.
if [[ -n "$CONDA_SH" && -f "$CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
  conda activate "$ENV_PREFIX"
  ENV_KIND="conda_sh"
elif [[ -x "$MICROMAMBA_BIN" ]]; then
  # POSIX activation script is used instead of a login shell so `set -e` stays intact.
  eval "$("$MICROMAMBA_BIN" shell posix hook --prefix "$ENV_PREFIX" 2>/dev/null)" || true
  micromamba activate "$ENV_PREFIX"
  ENV_KIND="micromamba"
elif [[ -f "$ENV_PREFIX/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "$ENV_PREFIX/bin/activate"
  ENV_KIND="venv_activate"
elif [[ -x "$ENV_PREFIX/bin/python" ]]; then
  export PATH="$ENV_PREFIX/bin:$PATH"
  ENV_KIND="path_prefix"
else
  echo "Unable to activate environment: $ENV_PREFIX" >&2
  exit 2
fi
export ULTIMATE_ENV_ACTIVATION="$ENV_KIND"
echo "[ultimate_run] env_prefix=$ENV_PREFIX activation=$ENV_KIND runtime_root=$RUNTIME_ROOT" >&2
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  export ULTIMATE_THREADS="${ULTIMATE_THREADS:-$SLURM_CPUS_PER_TASK}"
fi
cd "$PROJECT_ROOT"
if [[ -n "$PRODUCTION_APPROVAL" ]]; then
  PRODUCTION_APPROVAL="$(readlink -f "$PRODUCTION_APPROVAL")"
  PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}" python -m ultimate.cli run --config "$CONFIG_PATH" --production-approval "$PRODUCTION_APPROVAL"
else
  PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}" python -m ultimate.cli run --config "$CONFIG_PATH"
fi
