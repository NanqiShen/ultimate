#!/usr/bin/env bash
set -euo pipefail

MODE="local"
ROOT=""
PYTHON_BIN=""

usage() {
  cat <<'USAGE'
usage: scripts/check_dev_entrypoint.sh [--mode local|remote] [--root PATH] [--python PATH]

Checks the Ultimate development entrypoint without requiring PYTHONPATH.
Remote mode runs the same lightweight checks on the shared root through SSH
(override the client with ULTIMATE_SSH) instead of a site-specific wrapper.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2:?--mode requires local or remote}"
      shift 2
      ;;
    --root)
      ROOT="${2:?--root requires a path}"
      shift 2
      ;;
    --python)
      PYTHON_BIN="${2:?--python requires a path}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "$MODE" != "local" && "$MODE" != "remote" ]]; then
  echo "--mode must be local or remote" >&2
  exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_local_root="$(cd "$script_dir/.." && pwd)"
if [[ -z "$ROOT" ]]; then
  if [[ "$MODE" == "remote" ]]; then
    ROOT="${ULTIMATE_ROOT:?Specify --root or set ULTIMATE_ROOT for remote mode}"
  else
    ROOT="${ULTIMATE_ROOT:-$default_local_root}"
  fi
fi

# Managed environments live outside the code root on shared deployments, so a
# site points at them explicitly instead of relying on <root>/.conda/envs.
ENV_ROOT="${ULTIMATE_ENV_ROOT:-$(dirname "$ROOT")/ultimate_runtime/envs}"

run_local_checks() {
  local root="$1"
  local python_bin="${2:-python}"
  local tmp_dir
  tmp_dir="$(mktemp -d)"
  trap "rm -rf '$tmp_dir'" EXIT

  cd "$root"
  echo "[ultimate-dev-entrypoint] root=$root"
  echo "[ultimate-dev-entrypoint] python=$python_bin"
  "$python_bin" - <<'PY'
import sys
if sys.version_info < (3, 11):
    raise SystemExit(f"Ultimate requires Python >=3.11 for default dev checks; got {sys.version.split()[0]}")
PY
  "$python_bin" -m pytest -q \
    tests/test_cli.py::test_cli_prepare_intake \
    tests/test_analysis_levels.py \
    tests/test_packaging_stability.py \
    tests/test_handoff_check.py
  "$python_bin" -m pip install -e .
  local ultimate_bin
  ultimate_bin="$(dirname "$("$python_bin" -c 'import sys; print(sys.executable)')")/ultimate"
  if [[ ! -x "$ultimate_bin" ]]; then
    ultimate_bin="$(command -v ultimate)"
  fi
  "$ultimate_bin" --help >/dev/null
  "$ultimate_bin" handoff-check --root "$root" --output-dir "$tmp_dir/handoff_check" >/dev/null
  "$ultimate_bin" prepare-intake --root "$tmp_dir/root" --output-dir "$tmp_dir/intake" --refresh-audit >/dev/null
  "$python_bin" - "$tmp_dir" <<'PY'
import json
import sys
from pathlib import Path

tmp = Path(sys.argv[1])
handoff = json.loads((tmp / "handoff_check" / "handoff_check_manifest.json").read_text(encoding="utf-8"))
if handoff.get("status") != "ready":
    raise SystemExit(f"handoff-check not ready: {handoff}")
intake = json.loads((tmp / "intake" / "intake_package_manifest.json").read_text(encoding="utf-8"))
template_status = intake.get("template_status") if isinstance(intake.get("template_status"), dict) else {}
if template_status.get("status") != "ready":
    raise SystemExit(f"prepare-intake template lookup not ready: {intake}")
print("dev_entrypoint_status=ready")
PY
}

run_remote_checks() {
  local root="$1"
  local python_bin="$2"
  local ssh_cmd="${ULTIMATE_SSH:-ssh}"
  if [[ -z "$python_bin" ]]; then
    if [[ -x "$ENV_ROOT/ultimate-core/bin/python" ]]; then
      python_bin="$ENV_ROOT/ultimate-core/bin/python"
    else
      python_bin="$(command -v python3)"
    fi
  fi
  if ! command -v "${ssh_cmd%% *}" >/dev/null 2>&1; then
    echo "--mode remote requires an SSH client; set ULTIMATE_SSH or use --mode local" >&2
    exit 2
  fi
  # shellcheck disable=SC2086
  $ssh_cmd bash -s -- "$root" "$python_bin" "$ENV_ROOT" <<'REMOTE'
set -euo pipefail
ROOT="$1"
PYTHON_BIN="$2"
ENV_ROOT="$3"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
cd "$ROOT"
echo "[ultimate-dev-entrypoint] root=$ROOT"
echo "[ultimate-dev-entrypoint] python=$PYTHON_BIN"
"$PYTHON_BIN" -m pytest -q \
  tests/test_cli.py::test_cli_prepare_intake \
  tests/test_analysis_levels.py \
  tests/test_packaging_stability.py \
  tests/test_handoff_check.py
"$PYTHON_BIN" -m pip install -e .
ULTIMATE_BIN="$ENV_ROOT/ultimate-core/bin/ultimate"
if [[ ! -x "$ULTIMATE_BIN" ]]; then
  ULTIMATE_BIN="$(command -v ultimate)"
fi
"$ULTIMATE_BIN" --help >/dev/null
"$ULTIMATE_BIN" handoff-check --root "$ROOT" --output-dir "$TMP_DIR/handoff_check" >/dev/null
"$ULTIMATE_BIN" prepare-intake --root "$TMP_DIR/root" --output-dir "$TMP_DIR/intake" --refresh-audit >/dev/null
"$PYTHON_BIN" - "$TMP_DIR" <<'PY'
import json
import sys
from pathlib import Path

tmp = Path(sys.argv[1])
handoff = json.loads((tmp / "handoff_check" / "handoff_check_manifest.json").read_text(encoding="utf-8"))
if handoff.get("status") != "ready":
    raise SystemExit(f"handoff-check not ready: {handoff}")
intake = json.loads((tmp / "intake" / "intake_package_manifest.json").read_text(encoding="utf-8"))
template_status = intake.get("template_status") if isinstance(intake.get("template_status"), dict) else {}
if template_status.get("status") != "ready":
    raise SystemExit(f"prepare-intake template lookup not ready: {intake}")
print("dev_entrypoint_status=ready")
PY
REMOTE
}

if [[ "$MODE" == "remote" ]]; then
  run_remote_checks "$ROOT" "$PYTHON_BIN"
else
  run_local_checks "$ROOT" "${PYTHON_BIN:-python}"
fi
