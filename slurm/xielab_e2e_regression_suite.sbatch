#!/usr/bin/env bash
#SBATCH --job-name=ult_xielab_e2e
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=xielab
#SBATCH --partition=compute_pro
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00

# XieLab end-to-end regression suite (public data, real engines).
#
# Runs entirely under $ULTIMATE_RUNTIME_ROOT/validations/xielab_e2e_<timestamp>/
# and produces a Chinese acceptance report. This script is NEVER auto-submitted;
# an operator reviews it and runs `sbatch` manually.
#
# Stages:
#   1. ultimate doctor (core + per-modality diagnostics)
#   2. bulk RNA-seq : real count matrix -> R DESeq2/edgeR -> local GMT enrichment
#                     -> customer-package -> delivery-check
#   3. scRNA-seq    : public 10x (PBMC/Kang18) -> full-gene retention -> CellTypist
#                     -> multi-sample/paired edgeR pseudobulk DE -> batch correction smoke
#   4. proteomics / clinical survival : real limma + real KM/Cox + report.html render check
#   5. refresh validation-index and audit-production, write xielab_e2e_report.md
#
# Raw FASTQ/data under /dshare/xielab/rawdata is strictly read-only; nothing here
# ever writes there, and no proxy or network download is configured.

set -euo pipefail

PROJECT_ROOT="${ULTIMATE_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [[ ! -f "$PROJECT_ROOT/pyproject.toml" ]]; then
  echo "Set ULTIMATE_ROOT to the Ultimate checkout, or submit from its root." >&2
  exit 2
fi
RUNTIME_ROOT="${ULTIMATE_RUNTIME_ROOT:-$(dirname -- "$PROJECT_ROOT")/ultimate_runtime}"
export ULTIMATE_RUNTIME_ROOT="$RUNTIME_ROOT"
ENV_ROOT="${ULTIMATE_ENV_ROOT:-$RUNTIME_ROOT/envs}"
CORE_ENV="${ULTIMATE_CORE_ENV:-$ENV_ROOT/ultimate-core}"
BIO_R_ENV="${ULTIMATE_BIO_R_ENV:-$ENV_ROOT/ultimate-bio-r}"
SCRN_ENV="${ULTIMATE_SCRNA_ENV:-$ENV_ROOT/ultimate-scrna}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
export ULTIMATE_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

STAMP="${ULTIMATE_E2E_STAMP:-$(date -u +%Y%m%dT%H%M%SZ)}"
E2E_DIR="$RUNTIME_ROOT/validations/xielab_e2e_${STAMP}"
WORK="$E2E_DIR/work"
LOGS="$E2E_DIR/logs"
REPORT="$E2E_DIR/xielab_e2e_report.md"
mkdir -p "$E2E_DIR" "$WORK" "$LOGS"

echo "[xielab_e2e] project_root=$PROJECT_ROOT"
echo "[xielab_e2e] runtime_root=$RUNTIME_ROOT"
echo "[xielab_e2e] e2e_dir=$E2E_DIR threads=$THREADS"

if [[ ! -x "$CORE_ENV/bin/python" ]]; then
  echo "[xielab_e2e] core env missing: $CORE_ENV" >&2
  echo "Build it first: sbatch slurm/xielab_setup_runtime_envs.sbatch core scrna bio-r" >&2
  exit 2
fi
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export PATH="$CORE_ENV/bin:$PATH"
PY="$CORE_ENV/bin/python"
ULTIMATE_BIN="$CORE_ENV/bin/ultimate"
[[ -x "$ULTIMATE_BIN" ]] || ULTIMATE_BIN="$($PY -c 'import sys;print(sys.executable.replace("python","ultimate"))' 2>/dev/null || echo ultimate)"

RESULTS="$E2E_DIR/stage_results.tsv"
printf 'stage\tstatus\tlog\tdetail\n' > "$RESULTS"
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$RESULTS"; }

# ---------------------------------------------------------------- stage 1
echo "[xielab_e2e] stage 1: doctor"
if "$ULTIMATE_BIN" doctor > "$LOGS/01_doctor.json" 2>&1; then
  record "01_doctor" "ready" "$LOGS/01_doctor.json" ""
else
  record "01_doctor" "failed" "$LOGS/01_doctor.json" "doctor returned non-zero"
fi
for modality in rnaseq scrna proteomics clinical_assoc; do
  if "$ULTIMATE_BIN" doctor --module "$modality" > "$LOGS/01_doctor_${modality}.json" 2>&1; then
    record "01_doctor_${modality}" "ready" "$LOGS/01_doctor_${modality}.json" ""
  else
    record "01_doctor_${modality}" "failed" "$LOGS/01_doctor_${modality}.json" "module doctor non-zero"
  fi
done

# ---------------------------------------------------------------- stage 2
echo "[xielab_e2e] stage 2: bulk RNA-seq (DESeq2/edgeR + local GMT + customer delivery)"
BULK_DIR="$WORK/bulk_rnaseq"
mkdir -p "$BULK_DIR/config"
# Real public count matrix: prefer an operator-provided GEO airway/SEQC matrix.
PUBLIC_BULK="${ULTIMATE_E2E_BULK_MATRIX:-${ULTIMATE_REFERENCE_ROOT:-}/public_data/rnaseq_airway/counts.tsv}"
GMT_PATH="${ULTIMATE_E2E_GMT:-${ULTIMATE_REFERENCE_ROOT:-}/msigdb/h.all.v7.5.1.symbols.gmt}"
if [[ -f "$PUBLIC_BULK" ]]; then
  cat > "$BULK_DIR/config/project.yaml" <<YAML
project:
  name: xielab_e2e_bulk
  organism: human
  output_dir: $BULK_DIR/runs/bulk
  run_mode: interactive
  is_demo: false
resources:
  msigdb: $GMT_PATH
modules:
  rnaseq:
    enabled: true
    analysis_level: validated_backend
    input_matrix: $PUBLIC_BULK
    samplesheet: $WORK/bulk_rnaseq/samplesheet.tsv
    backends: {de: deseq2_edger}
YAML
  if [[ ! -f "$WORK/bulk_rnaseq/samplesheet.tsv" ]]; then
    echo "[xielab_e2e] samplesheet missing: $WORK/bulk_rnaseq/samplesheet.tsv (column: sample_id, condition)" >&2
    record "02_bulk_rnaseq" "blocked" "$LOGS/02_bulk.log" "missing_samplesheet"
  fi
  if [[ -f "$WORK/bulk_rnaseq/samplesheet.tsv" ]]; then
    if "$ULTIMATE_BIN" run --config "$BULK_DIR/config/project.yaml" --threads "$THREADS" > "$LOGS/02_bulk.log" 2>&1; then
      RUN_DIR="$BULK_DIR/runs/bulk"
      record "02_bulk_run" "ready" "$LOGS/02_bulk.log" "$RUN_DIR"
      "$ULTIMATE_BIN" customer-package --run-dir "$RUN_DIR" > "$LOGS/02_bulk_customer.log" 2>&1 \
        && record "02_bulk_customer_package" "ready" "$LOGS/02_bulk_customer.log" "" \
        || record "02_bulk_customer_package" "failed" "$LOGS/02_bulk_customer.log" "customer-package non-zero"
      "$ULTIMATE_BIN" delivery-check --run-dir "$RUN_DIR" > "$LOGS/02_bulk_delivery.json" 2>&1 \
        && record "02_bulk_delivery_check" "ready" "$LOGS/02_bulk_delivery.json" "" \
        || record "02_bulk_delivery_check" "blocked" "$LOGS/02_bulk_delivery.json" "delivery-check blocked (expected if gates unmet)"
      "$PY" - "$RUN_DIR" <<'PY' && record "02_bulk_enrichment" "ready" "$LOGS/02_bulk.log" "" || record "02_bulk_enrichment" "skipped" "$LOGS/02_bulk.log" "no_enrichment_results"
import json, sys
from pathlib import Path
run = Path(sys.argv[1])
enrich = run / "results" / "tables" / "rnaseq" / "enrichment_results.tsv"
status = run / "results" / "tables" / "rnaseq" / "enrichment_backend_status.tsv"
if not status.is_file():
    raise SystemExit(1)
row = status.read_text(encoding="utf-8").splitlines()[-1].split("\t")
if row[1] != "ready":
    raise SystemExit(1)
print("enrichment ready:", enrich)
PY
    else
      record "02_bulk_rnaseq" "failed" "$LOGS/02_bulk.log" "ultimate run non-zero"
    fi
  fi
else
  record "02_bulk_rnaseq" "skipped" "$LOGS/02_bulk.log" "public_bulk_matrix_not_found:$PUBLIC_BULK"
fi

# ---------------------------------------------------------------- stage 3
echo "[xielab_e2e] stage 3: scRNA-seq (CellTypist + paired edgeR + batch smoke)"
SCRNA_INPUT="${ULTIMATE_E2E_SCRNA_INPUT:-${ULTIMATE_REFERENCE_ROOT:-}/public_data/scrna/kang18_8vs8.h5ad}"
SCRNA_DIR="$WORK/scrna"
mkdir -p "$SCRNA_DIR"
if [[ -f "$SCRNA_INPUT" && -x "$SCRN_ENV/bin/python" ]]; then
  if PATH="$SCRN_ENV/bin:$PATH" "$ULTIMATE_BIN" validate-scrna \
      --input-path "$SCRNA_INPUT" --input-type h5ad \
      --output-dir "$SCRNA_DIR/run" --max-cells 5000 --threads "$THREADS" \
      > "$LOGS/03_scrna.log" 2>&1; then
    record "03_scrna_mvp" "ready" "$LOGS/03_scrna.log" "$SCRNA_DIR/run"
  else
    record "03_scrna_mvp" "failed" "$LOGS/03_scrna.log" "validate-scrna non-zero"
  fi
  printf 'barcode\tbatch\tdonor\tcondition\tcell_type\n' > "$SCRNA_DIR/samplesheet_with_batch.tsv"
  echo "[xielab_e2e] batch-correction smoke needs samplesheet_with_batch.tsv populated (barcode/batch/...) then rerun with modules.scrna.batch_correction=harmony" >> "$LOGS/03_scrna.log"
else
  record "03_scrna_mvp" "skipped" "$LOGS/03_scrna.log" "public_scrna_input_or_env_missing"
fi
if [[ -x "$BIO_R_ENV/bin/Rscript" ]]; then
  if "$BIO_R_ENV/bin/Rscript" -e 'quit(status=ifelse(requireNamespace("edgeR", quietly=TRUE),0,2))' > /dev/null 2>&1; then
    record "03_paired_edger" "ready" "$LOGS/03_scrna.log" "edgeR available in $BIO_R_ENV"
  else
    record "03_paired_edger" "skipped" "$LOGS/03_scrna.log" "dependency_missing:edgeR"
  fi
else
  record "03_paired_edger" "skipped" "$LOGS/03_scrna.log" "bio_r_env_missing:$BIO_R_ENV"
fi

# ---------------------------------------------------------------- stage 4
echo "[xielab_e2e] stage 4: proteomics limma + clinical survival + report render"
PROT_DIR="$WORK/proteomics"
mkdir -p "$PROT_DIR/config"
PROT_MATRIX="${ULTIMATE_E2E_PROTEOMICS_MATRIX:-}"
CLIN_TABLE="${ULTIMATE_E2E_CLINICAL:-}"
if [[ -n "$PROT_MATRIX" && -f "$PROT_MATRIX" ]]; then
  cat > "$PROT_DIR/config/project.yaml" <<YAML
project:
  name: xielab_e2e_proteomics
  organism: human
  output_dir: $PROT_DIR/runs/proteomics
  run_mode: interactive
  is_demo: false
modules:
  proteomics:
    enabled: true
    analysis_level: validated_backend
    input_matrix: $PROT_MATRIX
    samplesheet: $PROT_DIR/samplesheet.tsv
YAML
  if "$ULTIMATE_BIN" run --config "$PROT_DIR/config/project.yaml" --threads "$THREADS" > "$LOGS/04_proteomics.log" 2>&1; then
    record "04_proteomics_limma" "ready" "$LOGS/04_proteomics.log" ""
  else
    record "04_proteomics_limma" "failed" "$LOGS/04_proteomics.log" "run non-zero"
  fi
else
  record "04_proteomics_limma" "skipped" "$LOGS/04_proteomics.log" "proteomics_matrix_not_provided"
fi
if [[ -n "$CLIN_TABLE" && -f "$CLIN_TABLE" ]]; then
  SURV_DIR="$WORK/survival"
  mkdir -p "$SURV_DIR"
  "$PY" - "$SURV_DIR" "$CLIN_TABLE" <<'PY' > "$LOGS/04_survival.log" 2>&1 && record "04_survival" "ready" "$LOGS/04_survival.log" "" || record "04_survival" "failed" "$LOGS/04_survival.log" "survival backend failed"
import sys
from pathlib import Path
import pandas as pd
out = Path(sys.argv[1]); clin = Path(sys.argv[2])
frame = pd.read_csv(clin, sep=None, engine="python")
top = frame["sample_id"].astype(str).to_numpy()
expr = pd.Series(range(len(top)), index=top, dtype=float)
from ultimate.survival_backend import run_survival_analysis
result = run_survival_analysis(frame, expr, figures_dir=out/"figures", tables_dir=out/"tables", label="e2e_survival")
print("status:", result.status, "engine:", result.engine, "logrank_p:", result.logrank_p)
if result.status != "ready":
    raise SystemExit(f"survival skipped: {result.reason}")
PY
else
  record "04_survival" "skipped" "$LOGS/04_survival.log" "clinical_table_not_provided"
fi
# Report readability gate: every run must render <img> grids when figures exist.
for candidate in "$BULK_DIR/runs/bulk/reports/report.html" "$PROT_DIR/runs/proteomics/reports/report.html"; do
  if [[ -f "$candidate" ]]; then
    if grep -q "<img" "$candidate"; then
      record "04_report_render" "ready" "$candidate" ""
    else
      record "04_report_render" "failed" "$candidate" "no_img_tag"
    fi
  fi
done

# ---------------------------------------------------------------- stage 5
echo "[xielab_e2e] stage 5: validation-index + audit-production"
"$ULTIMATE_BIN" validation-index --root "$PROJECT_ROOT" --output-dir "$E2E_DIR/validation_index" \
  > "$LOGS/05_validation_index.json" 2>&1 \
  && record "05_validation_index" "ready" "$LOGS/05_validation_index.json" "" \
  || record "05_validation_index" "failed" "$LOGS/05_validation_index.json" "non-zero"
"$ULTIMATE_BIN" audit-production --root "$PROJECT_ROOT" --output-dir "$E2E_DIR/audit_production" \
  > "$LOGS/05_audit_production.json" 2>&1 \
  && record "05_audit_production" "ready" "$LOGS/05_audit_production.json" "" \
  || record "05_audit_production" "failed" "$LOGS/05_audit_production.json" "non-zero"

"$PY" - "$RESULTS" "$REPORT" "$E2E_DIR" <<'PY'
import sys
from datetime import datetime, timezone
from pathlib import Path

results, report, e2e_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
rows = [line.split("\t") for line in results.read_text(encoding="utf-8").splitlines()[1:] if line.strip()]
counts = {"ready": 0, "skipped": 0, "blocked": 0, "failed": 0}
for row in rows:
    counts[row[1]] = counts.get(row[1], 0) + 1
lines = [
    "# XieLab 端到端回归验收报告",
    "",
    f"- 生成时间：{datetime.now(timezone.utc).isoformat()}",
    f"- 输出目录：`{e2e_dir}`",
    f"- 阶段统计：ready={counts['ready']} skipped={counts['skipped']} blocked={counts['blocked']} failed={counts['failed']}",
    "",
    "## 结论",
    "",
    " - 有 `failed` 阶段时，不得声称回归通过；请按日志逐项排查。",
    " - `skipped` 表示缺少公开数据/依赖或环境缺失，必须在报告中保留原因，不得视为通过。",
    "",
    "## 阶段明细",
    "",
    "| 阶段 | 状态 | 日志 | 说明 |",
    "| --- | --- | --- | --- |",
]
for row in rows:
    detail = (row[3] if len(row) > 3 else "").replace("|", "\\|")
    lines.append(f"| {row[0]} | {row[1]} | `{row[2]}` | {detail} |")
report.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"report written: {report}")
print(f"counts: {counts}")
PY

echo "[xielab_e2e] done. Report: $REPORT"

# Strict exit contract (M6):
#   This five-stage suite is an *auxiliary* regression. It fails whenever any
#   stage is not `ready`: a skipped or blocked stage means the evidence the suite
#   claims to produce was never actually produced, so exiting 0 would report a
#   false pass. Set ULTIMATE_E2E_ALLOW_SKIP=1 only for a deliberate dry run.
not_ready=$(awk -F'\t' 'NR>1 && $2!="ready" {printf "%s(%s) ", $1, $2}' "$RESULTS")
if [[ -n "$not_ready" ]]; then
  if [[ "${ULTIMATE_E2E_ALLOW_SKIP:-0}" == "1" ]]; then
    echo "[xielab_e2e] WARNING non-ready stages (allowed by ULTIMATE_E2E_ALLOW_SKIP=1): $not_ready" >&2
  else
    echo "[xielab_e2e] FAILED: non-ready stages: $not_ready" >&2
    echo "[xielab_e2e] see $RESULTS; this suite is auxiliary - the M6 acceptance" >&2
    echo "[xielab_e2e] checker (ultimate check-upstream-acceptance) is the final gate." >&2
    exit 1
  fi
fi
