# nf-core/sopa Handoff Template

用途：参考 nf-core/sopa / SOPA 生态处理 Xenium、CosMX、MERSCOPE、Visium HD 等空间组学，输出 SpatialData/SOPA project，再交回 Ultimate `spatial`。

本模板是 handoff，不表示 Ultimate 自动完成所有厂商格式解析。

## 命令计划

```bash
set -euo pipefail

export NXF_HOME=${ULTIMATE_ROOT}/.nextflow
export NXF_WORK=${ULTIMATE_ROOT}/jobs/<job_id>/work/nfcore_sopa

nextflow run nf-core/sopa \
  -profile apptainer,slurm \
  -c ${ULTIMATE_ROOT}/jobs/<job_id>/config/nfcore_sopa_nextflow.config \
  -params-file ${ULTIMATE_ROOT}/jobs/<job_id>/config/nfcore_sopa_params.yaml \
  -work-dir "$NXF_WORK" \
  -resume
```

## Ultimate handback

- SpatialData Zarr 或 SOPA project path
- cell/spot metadata
- segmentation/coordinate tables
- QC summary and figures

Xenium/CosMX/MERSCOPE 不能强行走 Visium 逻辑；图像和 segmentation 质量必须单独写限制。
