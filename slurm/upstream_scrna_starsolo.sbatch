#!/usr/bin/env bash
#SBATCH --job-name=ult_up_scrna
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=xielab
#SBATCH --partition=compute_pro
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00

# Standard open-source scRNA-seq upstream: STAR --soloType CB_UMI_Simple ->
# standard 10x `filtered/` directory + import_config.yaml (Ultimate 10x_mtx import).
#
# Generated/reviewed by `ultimate prepare-upstream`; raw FASTQ stays read-only and
# no network access is required. Threads are bound to $SLURM_CPUS_PER_TASK.
#
# Expected jobspec variables:
#   SAMPLESHEET  MANIFEST  FASTQ_DIR  INDEX  OUT_DIR  THREADS
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

if [[ -x "$UPSTREAM_ENV/bin/STAR" ]]; then
  export PATH="$UPSTREAM_ENV/bin:$PATH"
else
  echo "[upstream_scrna] upstream environment not found at $UPSTREAM_ENV" >&2
  echo "Build it with: sbatch slurm/xielab_setup_upstream_env.sbatch" >&2
  exit 2
fi

: "${SAMPLESHEET:?SAMPLESHEET must be set (FASTQ samplesheet.tsv)}"
: "${INDEX:?INDEX must be set (STAR genome index directory)}"
: "${OUT_DIR:?OUT_DIR must be set}"
FASTQ_DIR="${FASTQ_DIR:-$(dirname -- "$SAMPLESHEET")}"
mkdir -p "$OUT_DIR"/{logs,star,tables}

if [[ "$(head -1 "$SAMPLESHEET")" != $'sample_id\tfastq_1\tfastq_2\tstrandedness'* ]]; then
  echo "[upstream_scrna] SAMPLESHEET is not the normalised table: $SAMPLESHEET" >&2
  echo "Regenerate the package with: ultimate prepare-upstream ..." >&2
  exit 2
fi
echo "[upstream_scrna] samplesheet=$SAMPLESHEET (normalised)"
echo "[upstream_scrna] index=$INDEX"
echo "[upstream_scrna] out_dir=$OUT_DIR threads=$THREADS"

for required in "$SAMPLESHEET"; do
  if [[ ! -r "$required" ]]; then
    echo "[upstream_scrna] cannot read: $required" >&2
    exit 2
  fi
done
if [[ ! -d "$INDEX" ]]; then
  echo "[upstream_scrna] STAR index directory missing: $INDEX" >&2
  exit 2
fi
if [[ ! -f "$INDEX/SAindex" && ! -f "$INDEX/Genome" ]]; then
  echo "[upstream_scrna] INDEX does not look like a STAR index (missing SAindex/Genome): $INDEX" >&2
  exit 2
fi

# STAR filters must match the whitelist shipped with the configured chemistry.
WL="${WHITELIST:-}"
if [[ -z "$WL" ]]; then
  for candidate in \
    "$UPSTREAM_ENV/share/star/3M-february-2018.txt" \
    "$UPSTREAM_ENV/lib/python"*/site-packages/whitelist/3M-february-2018.txt
  do
    if [[ -f "$candidate" ]]; then WL="$candidate"; break; fi
  done
fi
if [[ -z "$WL" || ! -r "$WL" ]]; then
  echo "[upstream_scrna] 10x barcode whitelist not found; set WHITELIST=<path to 3M-february-2018.txt>" >&2
  exit 2
fi

: > "$OUT_DIR/tables/import_configs.txt"
while IFS=$'\t' read -r sample r1 r2 _strandedness; do
  [[ -z "$sample" ]] && continue
  case "$r1" in /*) fq1="$r1" ;; *) fq1="$FASTQ_DIR/$r1" ;; esac
  if [[ ! -r "$fq1" ]]; then
    echo "[upstream_scrna] R1 not readable: $fq1" >&2
    exit 2
  fi
  case "$r2" in
    /*) fq2="$r2" ;;
    "") fq2="" ;;
    *) fq2="$FASTQ_DIR/$r2" ;;
  esac
  if [[ -n "$fq2" && ! -r "$fq2" ]]; then
    echo "[upstream_scrna] R2 not readable: $fq2" >&2
    exit 2
  fi
  sample_out="$OUT_DIR/star/${sample}"
  mkdir -p "$sample_out"
  reads=("$fq1")
  [[ -n "$fq2" ]] && reads+=("$fq2")
  echo "[upstream_scrna] STARsolo $sample"
  STAR --runThreadN "$THREADS" \
       --genomeDir "$INDEX" \
       --readFilesIn "${reads[@]}" \
       --readFilesCommand zcat \
       --soloType CB_UMI_Simple \
       --soloCBwhitelist "$WL" \
       --soloFeatures Gene \
       --soloCellFilter EmptyDrops_CR \
       --outSAMtype None \
       --outFileNamePrefix "$sample_out/" \
       > "$OUT_DIR/logs/${sample}.star.log" 2>&1

  # Standard 10x layout so Ultimate can ingest with input_type=10x_mtx.
  filtered="$sample_out/filtered"
  mkdir -p "$filtered"
  gene_dir="$sample_out/Solo.out/Gene/filtered"
  if [[ ! -f "$gene_dir/matrix.mtx" ]]; then
    echo "[upstream_scrna] STARsolo filtered matrix missing: $gene_dir/matrix.mtx" >&2
    exit 2
  fi
  cp "$gene_dir/matrix.mtx" "$filtered/matrix.mtx"
  if [[ -f "$gene_dir/barcodes.tsv" ]]; then
    cp "$gene_dir/barcodes.tsv" "$filtered/barcodes.tsv"
  fi
  if [[ -f "$gene_dir/features.tsv" ]]; then
    cp "$gene_dir/features.tsv" "$filtered/features.tsv"
  elif [[ -f "$sample_out/Solo.out/Gene/features.tsv" ]]; then
    cp "$sample_out/Solo.out/Gene/features.tsv" "$filtered/features.tsv"
  fi
  gzip -f "$filtered/matrix.mtx"
  [[ -f "$filtered/barcodes.tsv" ]] && gzip -f "$filtered/barcodes.tsv"
  [[ -f "$filtered/features.tsv" ]] && gzip -f "$filtered/features.tsv"

  cat > "$OUT_DIR/import_config_${sample}.yaml" <<YAML
module: scrna
input_type: 10x_mtx
input_path: $filtered
sample_id: $sample
upstream_route: open_source_starsolo
notes: STARsolo CB_UMI_Simple output reshaped into the standard 10x filtered/ layout.
YAML
  echo -e "${sample}\t$OUT_DIR/import_config_${sample}.yaml" >> "$OUT_DIR/tables/import_configs.txt"
done < <(tail -n +2 "$SAMPLESHEET")

python3 - "$OUT_DIR" <<'PY'
import json, sys
from pathlib import Path

out_dir = Path(sys.argv[1])
configs = [line.split("\t")[1] for line in (out_dir / "tables" / "import_configs.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
if not configs:
    raise SystemExit("no STARsolo import configs produced")
(out_dir / "import_config.yaml").write_text(
    "module: scrna\n"
    "input_type: 10x_mtx_multi_sample\n"
    "samples:\n" + "".join(f"  - config: {path}\n" for path in configs) +
    "upstream_route: open_source_starsolo\n",
    encoding="utf-8",
)
(out_dir / "raw_qc_manifest.json").write_text(
    json.dumps({
        "status": "ready",
        "route": "STARsolo_CB_UMI_Simple",
        "n_samples": len(configs),
        "sample_configs": configs,
    }, indent=2),
    encoding="utf-8",
)
print(f"[upstream_scrna] wrote {len(configs)} import configs")
PY

echo "[upstream_scrna] done: $OUT_DIR/import_config.yaml"
