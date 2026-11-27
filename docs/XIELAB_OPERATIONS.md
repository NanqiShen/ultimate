# XieLab 集群部署与运维手册

本手册描述 Ultimate 在 XieLab HPC（`xielab`）上的唯一部署形态、环境矩阵、Slurm 参数与存储维护规范。

## 1. 单源目录（唯一真相源）

| 用途 | 路径 | 说明 |
| --- | --- | --- |
| 唯一代码根 | `/dshare/xielab/analysisdata/TFC/shennq/ultimate` | 只放代码与模板，**不放假数据、不放环境** |
| 唯一运行时根 | `/dshare/xielab/analysisdata/TFC/shennq/ultimate_runtime` | `envs/`、`micromamba/`、`micromamba-root/`、`cache/`、`logs/`、`manifests/`、`validations/` |
| 原始数据根 | `/dshare/xielab/rawdata/TFC/shennq` | **绝对只读**（`/dshare/rawdata` 全树同理），任何脚本不得写入或拷贝 |

规则：
- 代码根与运行时根是**同级目录**，`paths.runtime_root()` 会自动解析 `<code>/../ultimate_runtime`。
- 禁止把项目或环境建在 `/dshare/dhome`（家目录配额小且不被备份）。
- 登录节点只做轻量检查（`ultimate doctor`、`handoff-check`、`preflight`），重计算一律 `sbatch`。

## 2. 环境矩阵（4 套）

所有环境由 `$ULTIMATE_RUNTIME_ROOT/micromamba/bin/micromamba` 创建，落在 `$ULTIMATE_RUNTIME_ROOT/envs/`：

| 环境名 | 定义文件 | 用途 | 关键内容 |
| --- | --- | --- | --- |
| `ultimate-core` | `envs/xielab.core.yml` | 主流程 / CLI / 轻量统计 | python 3.11、numpy/pandas/scipy/matplotlib、lifelines、r-base、limma/DESeq2/edgeR、r-survival |
| `ultimate-scrna` | `envs/xielab.scrna.py.yml` | 单细胞 | scanpy、anndata、harmonypy、celltypist、scrublet、leidenalg、squidpy、gseapy、decoupler、liana |
| `ultimate-bio-r` | `envs/xielab.bio.r.yml` | R 生物信息后端 | limma、DESeq2、edgeR、clusterProfiler、fgsea、WGCNA、survival、survminer |
| `ultimate-upstream` | `envs/xielab.upstream.yml` | 开源上游定量 | fastqc、multiqc、fastp、salmon、star、subread、samtools、alevin-fry |

### 环境解析顺序

`src/ultimate/paths.py` 与 `src/ultimate/env_router.py` 共同决定运行时位置：

1. `ULTIMATE_ENV_ROOT`（显式覆盖）→ 否则 `$ULTIMATE_RUNTIME_ROOT/envs`
2. `resolve_rscript(module)`：`module_cfg.rscript` → `ULTIMATE_R_ENV_<MODULE>` / `ULTIMATE_RSCRIPT` → `environment_root()/<env>/bin/Rscript` → `PATH`
3. `resolve_packaged_asset(*parts)`：`ULTIMATE_ASSET_ROOT` → `ULTIMATE_ROOT` → 源码检出 → `sys.prefix/share/ultimate`

**显式配置了无效路径会直接报错，不会静默换用别的解释器。**

### 建环境（人工提交，不自动执行）

```bash
# 3 套运行时环境（core / scrna / bio-r）
ULTIMATE_ROOT=/dshare/xielab/analysisdata/TFC/shennq/ultimate \
  sbatch --account=xielab --partition=compute_pro --qos=normal \
  slurm/xielab_setup_runtime_envs.sbatch core scrna bio-r

# 开源上游环境（fastp/salmon/STAR/...）
ULTIMATE_ROOT=/dshare/xielab/analysisdata/TFC/shennq/ultimate \
  sbatch --account=xielab --partition=compute_pro --qos=normal \
  slurm/xielab_setup_upstream_env.sbatch
```

## 3. Slurm 参数（固定）

| 参数 | 值 |
| --- | --- |
| `account` | `xielab` |
| `partition` | `compute_pro` |
| `qos` | `normal` |
| sbatch 二进制 | 原生 `/usr/bin/sbatch`（`ULTIMATE_SBATCH` 可覆盖） |

- 集群**没有** `hpc-sbatch` / `hpc-run`；相关历史脚本已全部改为原生 `sbatch`。
- nextflow 模板的队列读 `System.getenv('ULTIMATE_SLURM_PARTITION') ?: 'compute_pro'`，不再硬编码 `comp`/`cpu`。
- 线程一律绑定 `$SLURM_CPUS_PER_TASK`：`ultimate run --threads`、生成脚本的 `#SBATCH --cpus-per-task`、`THREADS` 变量三者一致。

## 4. 开源上游路线（无 nextflow/apptainer）

`compute_pro` 未预装 nextflow / apptainer / cellranger，改走开源路线：

```bash
# 1) 生成可复核作业包（不会自动提交）
ultimate prepare-upstream --route salmon \
  --samplesheet <fastq_samplesheet.tsv> \
  --index <salmon_index_dir> \
  --output-dir jobs/<job_id>/config/upstream

# 2) 人工审阅后提交
sbatch jobs/<job_id>/config/upstream/run_upstream.sbatch
```

- `--route salmon`：`fastp` → `salmon quant` → `counts_raw.tsv` + `import_config.yaml` + `raw_qc_manifest.json`
- `--route starsolo`：`STAR --soloType CB_UMI_Simple` → 标准 10x `filtered/` + `import_config.yaml`
- 只读校验：FASTQ 必须存在且可读、扩展名为 `.fastq(.gz)/.fq(.gz)`；索引目录必须含 salmon `info.json` 或 STAR `SAindex`
- 缺 samplesheet / 缺 FASTQ / 缺索引 → 明确 `blocked` 并给出原因，不生成脚本

## 5. 存储维护规范（500 GB 预算）

运行时根可用空间按 500 GB 规划，超限前必须清理：

| 目录 | 建议上限 | 清理策略 |
| --- | --- | --- |
| `envs/` | 120 GB | 只保留 4 套在用环境；旧版环境确认无引用后 `micromamba env remove -p` |
| `cache/`（含 `public_data/`、`conda_pkgs/`、`nextflow/`） | 200 GB | 按访问时间清理；`conda_pkgs` 可安全删除后重下 |
| `validations/` | 120 GB | 保留最近 3 次 `xielab_e2e_*`；其余归档为报告后删除中间产物 |
| `logs/` | 20 GB | 保留 90 天，按 job id 轮转 |
| `manifests/` | 5 GB | 长期保留（体积小，是审计依据） |

日常检查：

```bash
ultimate audit-tools --root /dshare/xielab/analysisdata/TFC/shennq/ultimate   # 含存储估算
ultimate prune-tools --root ... --budget-gb 500                               # 生成清理计划（默认只出计划）
```

## 6. 端到端验收

```bash
ULTIMATE_ROOT=/dshare/xielab/analysisdata/TFC/shennq/ultimate \
  sbatch --account=xielab --partition=compute_pro --qos=normal \
  slurm/xielab_e2e_regression_suite.sbatch
```

套件在 `$ULTIMATE_RUNTIME_ROOT/validations/xielab_e2e_<timestamp>/` 下串联 doctor → bulk RNA-seq → scRNA → 蛋白/生存 → validation-index/audit-production，并输出中文报告 `xielab_e2e_report.md`。

判定规则（M6 起为严格契约）：
- 任一阶段不是 `ready`（含 `failed` / `blocked` / `skipped`）→ 套件**退出非 0**。
  仅在刻意演练时可用 `ULTIMATE_E2E_ALLOW_SKIP=1` 临时放行。
- 该五阶段套件只是**辅助回归**，M6 的最终判据是独立严格验收器：

```bash
ultimate check-upstream-acceptance \
  --package-dir <prepare-upstream 生成的包目录> \
  --downstream-manifest <下游 run_manifest.json>
```

验收器逐项检查：作业包产物、`upstream_out/` 实际矩阵路径、样本列与规范化样本表
一致、特征 ID 唯一、数值有限非负、非零信号、STARsolo 的 `filtered/` 目录维度与
feature/barcode 数相符、以及**下游确实消费了本次产物**。任一必需检查缺失、
跳过、阻塞或失败，整体退出非 0，并写出 `M6_ACCEPTANCE.json`。

### 计数语义（不得夸大）

- `--route salmon` 的 `counts_raw.tsv` 是 **salmon 转录本级 `NumReads` 四舍五入后的整数**，
  **不是基因级原始整数计数**；原始 `quant/<sample>/quant.sf` 必须保留，做基因级差异
  分析前需先做转录本→基因汇总。
- `--route starsolo` 输出的是 **STARsolo UMI 计数**（非负整数），格式为标准 10x
  `filtered/matrix.mtx(.gz)`。
