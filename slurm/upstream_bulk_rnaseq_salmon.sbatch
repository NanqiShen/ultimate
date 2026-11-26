#!/usr/bin/env bash
#SBATCH --job-name=ult_up_bulk
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=xielab
#SBATCH --partition=compute_pro
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# Standard open-source bulk RNA-seq upstream: fastp QC -> salmon quant ->
# Ultimate-importable counts_raw.tsv + import_config.yaml + raw_qc_manifest.json.
#
# This template is generated/reviewed by `ultimate prepare-upstream`; it never
# modifies raw FASTQ (raw inputs are strictly read-only) and never contacts the
# network. Threads are bound to $SLURM_CPUS_PER_TASK.
#
# Expected jobspec variables (written into run_upstream.sbatch by prepare-upstream):
#   SAMPLESHEET  MANIFEST  FASTQ_DIR  INDEX  OUT_DIR  THREADS  REFERENCE_LABEL
#
# Usage (after reviewing the generated copy):
#   sbatch jobs/<job_id>/config/run_upstream.sbatch

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
RUNTIME_ROOT="${ULTIMATE_RUNTIME_ROOT:-$(dirname -- "$PROJECT_ROOT")/ultimate_runtime}"
ENV_ROOT="${ULTIMATE_ENV_ROOT:-$RUNTIME_ROOT/envs}"
UPSTREAM_ENV="${ULTIMATE_UPSTREAM_ENV:-$ENV_ROOT/ultimate-upstream}"

THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-8}}"
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

if [[ -x "$UPSTREAM_ENV/bin/fastp" ]]; then
  export PATH="$UPSTREAM_ENV/bin:$PATH"
else
  echo "[upstream_bulk] upstream environment not found at $UPSTREAM_ENV" >&2
  echo "Build it with: sbatch slurm/xielab_setup_upstream_env.sbatch" >&2
  exit 2
fi

: "${SAMPLESHEET:?SAMPLESHEET must be set (FASTQ samplesheet.tsv)}"
: "${INDEX:?INDEX must be set (salmon index directory)}"
: "${OUT_DIR:?OUT_DIR must be set}"
FASTQ_DIR="${FASTQ_DIR:-$(dirname -- "$SAMPLESHEET")}"
mkdir -p "$OUT_DIR"/{fastp,quant,logs,tables}

echo "[upstream_bulk] samplesheet=$SAMPLESHEET"
echo "[upstream_bulk] index=$INDEX"
echo "[upstream_bulk] out_dir=$OUT_DIR threads=$THREADS"

# Raw FASTQ is read-only: verify we can read but never write next to it.
if [[ ! -r "$SAMPLESHEET" ]]; then
  echo "[upstream_bulk] cannot read samplesheet: $SAMPLESHEET" >&2
  exit 2
fi
if [[ ! -d "$INDEX" ]]; then
  echo "[upstream_bulk] salmon index directory missing: $INDEX" >&2
  exit 2
fi
if [[ ! -f "$INDEX/info.json" && ! -f "$INDEX/versionInfo.json" ]]; then
  echo "[upstream_bulk] INDEX does not look like a salmon index (missing info.json): $INDEX" >&2
  exit 2
fi

# SAMPLESHEET is the normalised table produced by `ultimate prepare-upstream`
# (columns: sample_id / fastq_1 / fastq_2 / strandedness). The validator and this
# job therefore consume the identical file; no re-parsing of raw column order.
if [[ "$(head -1 "$SAMPLESHEET")" != $'sample_id\tfastq_1\tfastq_2\tstrandedness'* ]]; then
  echo "[upstream_bulk] SAMPLESHEET is not the normalised table: $SAMPLESHEET" >&2
  echo "Regenerate the package with: ultimate prepare-upstream ..." >&2
  exit 2
fi
cp "$SAMPLESHEET" "$OUT_DIR/tables/samples_parsed.tsv"
if [[ "$(wc -l < "$OUT_DIR/tables/samples_parsed.tsv")" -lt 2 ]]; then
  echo "[upstream_bulk] samplesheet has no data rows: $SAMPLESHEET" >&2
  exit 2
fi

: > "$OUT_DIR/tables/quant_files.txt"
while IFS=$'\t' read -r sample r1 r2 strandedness; do
  [[ -z "$sample" ]] && continue
  case "$r1" in /*) fq1="$r1" ;; *) fq1="$FASTQ_DIR/$r1" ;; esac
  if [[ ! -r "$fq1" ]]; then
    echo "[upstream_bulk] R1 not readable: $fq1" >&2
    exit 2
  fi
  fastp_args=(--in1 "$fq1" --out1 "$OUT_DIR/fastp/${sample}_R1.trim.fq.gz"
              --json "$OUT_DIR/fastp/${sample}.fastp.json" --html "$OUT_DIR/fastp/${sample}.fastp.html"
              --thread "$THREADS" --report_title "$sample")
  if [[ -n "${r2:-}" && "$r2" != "-" ]]; then
    case "$r2" in /*) fq2="$r2" ;; *) fq2="$FASTQ_DIR/$r2" ;; esac
    [[ -r "$fq2" ]] || { echo "[upstream_bulk] R2 not readable: $fq2" >&2; exit 2; }
    fastp_args+=(--in2 "$fq2" --out2 "$OUT_DIR/fastp/${sample}_R2.trim.fq.gz")
  fi
  echo "[upstream_bulk] fastp $sample"
  fastp "${fastp_args[@]}" 2>> "$OUT_DIR/logs/${sample}.fastp.log"

  salmon_args=(-p "$THREADS" --libType "A" -i "$INDEX" -1 "$OUT_DIR/fastp/${sample}_R1.trim.fq.gz")
  if [[ -f "$OUT_DIR/fastp/${sample}_R2.trim.fq.gz" ]]; then
    salmon_args+=(-2 "$OUT_DIR/fastp/${sample}_R2.trim.fq.gz")
  fi
  echo "[upstream_bulk] salmon quant $sample"
  salmon quant "${salmon_args[@]}" -o "$OUT_DIR/quant/${sample}" \
    > "$OUT_DIR/logs/${sample}.salmon.log" 2>&1
  echo -e "${sample}\t$OUT_DIR/quant/${sample}/quant.sf" >> "$OUT_DIR/tables/quant_files.txt"
done < "$OUT_DIR/tables/samples_parsed.tsv"

python3 - "$OUT_DIR" <<'PY'
import sys
from pathlib import Path
import pandas as pd

out_dir = Path(sys.argv[1])
quant_list = out_dir / "tables" / "quant_files.txt"
frames = []
for line in quant_list.read_text(encoding="utf-8").splitlines():
    if not line.strip():
        continue
    sample, path = line.split("\t")
    frame = pd.read_csv(path, sep="\t")[["Name", "NumReads"]].rename(columns={"NumReads": sample})
    frames.append(frame.set_index("Name"))
if not frames:
    raise SystemExit("no salmon quant.sf produced")
matrix = pd.concat(frames, axis=1).fillna(0.0)
matrix.index.name = "feature_id"
matrix.round().astype("int64").to_csv(out_dir / "counts_raw.tsv", sep="\t")
(out_dir / "import_config.yaml").write_text(
    "module: rnaseq\n"
    "input_matrix: " + str(out_dir / "counts_raw.tsv") + "\n"
    "matrix_type: raw_counts_from_salmon\n"
    "upstream_route: open_source_salmon\n"
    "notes: salmon transcript-level NumReads rounded to integers; map to gene level before formal DE.\n",
    encoding="utf-8",
)
import json
(out_dir / "raw_qc_manifest.json").write_text(
    json.dumps({
        "status": "ready",
        "route": "fastp_salmon",
        "n_samples": int(matrix.shape[1]),
        "n_features": int(matrix.shape[0]),
        "counts_matrix": str(out_dir / "counts_raw.tsv"),
        "threads": "see SLURM_CPUS_PER_TASK in logs",
    }, indent=2),
    encoding="utf-8",
)
print(f"[upstream_bulk] wrote counts_raw.tsv {matrix.shape}")
PY

echo "[upstream_bulk] done: $OUT_DIR/counts_raw.tsv"
