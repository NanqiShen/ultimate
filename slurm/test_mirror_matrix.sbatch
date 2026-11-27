#!/usr/bin/env bash
#SBATCH --job-name=ult_mirror_mat
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" || ! -d "$PROJECT_ROOT/src/ultimate" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
export ULTIMATE_ROOT="$PROJECT_ROOT"
mkdir -p "$PROJECT_ROOT/logs" "$PROJECT_ROOT/audits/mirror_test"

python - <<'PY'
import json
import socket
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

mirrors = {
    "tuna": "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/repodata.json.zst",
    "bfsu": "https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/linux-64/repodata.json.zst",
    "sjtu": "https://mirror.sjtu.edu.cn/anaconda/cloud/conda-forge/linux-64/repodata.json.zst",
    "aliyun": "https://mirrors.aliyun.com/anaconda/cloud/conda-forge/linux-64/repodata.json.zst",
}
rows = []
for name, url in mirrors.items():
    host = urllib.request.urlparse(url).hostname
    row = {"name": name, "url": url, "host": host}
    try:
        row["addresses"] = sorted({item[4][0] for item in socket.getaddrinfo(host, 443)})
    except Exception as exc:
        row["dns_error"] = f"{type(exc).__name__}: {exc}"
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=20) as response:
            row["http_ok"] = True
            row["http_status"] = response.status
    except Exception as exc:
        row["http_ok"] = False
        row["http_error"] = f"{type(exc).__name__}: {exc}"
    rows.append(row)

out = Path(str(__import__("pathlib").Path(__import__("os").environ["ULTIMATE_ROOT"]) / "audits/mirror_test/mirror_matrix.json"))
out.write_text(json.dumps({"generated_at": datetime.now(timezone.utc).isoformat(), "results": rows}, indent=2, ensure_ascii=False), encoding="utf-8")
print(out.read_text(encoding="utf-8"))
PY

