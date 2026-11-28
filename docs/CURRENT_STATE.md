# Ultimate 当前状态基线

更新时间：2026-11-28（M1–M6 全线验收通过基线）

## 服务器基线

- 正式代码根：`/dshare/xielab/analysisdata/TFC/shennq/ultimate`（当前合流同步至 `NanqiShen/ultimate:main`，保留集群本地未跟踪 `AGENTS.md`）
- 运行时根：`/dshare/xielab/analysisdata/TFC/shennq/ultimate_runtime`（当前总占用约 `36G` / `500G` 预算）
- 四套 `micromamba` 环境全部 `ready`：
  - `ultimate-upstream`（`213484`：`fastp 1.4.0`、`salmon 2.12.0`、`STAR 2.7.10b`、`samtools 1.25`、`fastqc 0.13.0`、`multiqc 1.35`、`featureCounts`、`alevin-fry` + `3M-february-2018.txt`）
  - `ultimate-core`（`213485`：Python `3.11.16` + `ultimate-bioinfo 0.1.1rc1` + `lifelines`/`scipy`）
  - `ultimate-bio-r`（`213485`：`limma`/`DESeq2`/`edgeR`/`WGCNA`/`clusterProfiler`/`fgsea`/`survival`/`survminer`）
  - `ultimate-scrna`（`213494`：`scanpy 1.11.5`、`anndata 0.12.19`、`harmonypy 2.1.0`、`squidpy 1.8.2`、`gseapy 1.3.2`、`celltypist 1.7.1`、`decoupler 2.2.0`、`liana 1.10.0`）
  - 集群运维提示：计算节点 `node68` 存在外网 DNS 解析故障，建环境等需要拉取远端包时应加 `--exclude=node68`
- GitHub Actions CI：已跑通 `hardcoded path scan`、`pytest (py3.11/py3.12)`、`wheel install smoke (py3.11/py3.12)`
- 已完成 validation guard 标准化：`${ULTIMATE_ROOT}/audits/validation_guard_latest/validation_guard_check.tsv`
- 验证总索引：`${ULTIMATE_ROOT}/reports/validation_index/validation_index.tsv`
- 生产审计：`${ULTIMATE_ROOT}/audits/production_readiness_latest/production_audit.json`
- 肿瘤单细胞专项验证：`${ULTIMATE_ROOT}/validations/slurm_tumor_sc_maynard_raw_counts`

## 已通过状态

- 22 个模块标准外壳均为 `ready`。
- 严格 guard 后 production audit 当前为 `ready_basic=22`，无模块级 `partial`。
- `genotype_demux` 已用 vireo/cellSNP 公开示例矩阵完成 validated_backend；`hto_demux` 已用 Seurat 12-HTO 公开示例矩阵完成 validated_backend；`perturb_seq` 已用 Adamson 2016 public Perturb-seq h5ad fixture 完成 validated_backend。
- validation guard 当前无 non-ready 项；validation-index 当前汇总 34 个验证 run。
- final acceptance 当前通过 14 项，上一轮由 Perturb-seq synthetic/demo evidence 造成的 partial 已关闭。
- `tumor_sc` 已完成 NSCLC raw-count Slurm 专项验证，`analysis_level=validated_backend`，`delivery_allowed=false`。
- `slurm_tumor_sc_maynard_raw_counts` 基于 Maynard raw h5ad 抽样 3000 细胞，输出 16 个表、4 张图、1 个 h5ad 对象，并纳入 validation index 和 production audit。
- inferCNV、CopyKAT、Seurat 在 `ultimate-scrna-r` 中可用；Maynard 输入通过 raw integer count gate，CopyKAT/inferCNV 完整后端当前被人工放行 gate 阻断，阻断原因写入 `backend_attempts.tsv`。
- inferCNV gene order 已通过 Lambrechts 参考 h5ad 映射，29634 个基因中 22138 个完成染色体坐标映射。
- 单细胞能力矩阵仍保留授权工具边界：Cell Ranger、Space Ranger、CIBERSORT 等只做路径检测或 adapter，不默认声明为 fully automatic。

## 当前边界

- `validated_backend` 只代表公开数据或内部验证数据跑通，不等于客户正式交付。
- `production_backend` 必须经过 approval gate。
- demo/synthetic/smoke 结果不得作为正式交付证据。
- 高级工具若只完成输入准备或可用性检测，必须写成 `handoff`、`optional backend` 或 `adapter`。

## M1–M5 合并后能力状态（2026-11）

- **谱系与路径层（M1）**：XieLab TFC 7 个提交已合入统一基线；`paths.py` 同时提供
  `project_root/resource_root/external_root` 与 `runtime_root/environment_root/cache_root`，
  `default_platform_root`/`default_environment_root` 保留为别名。代码根与 `ultimate_runtime`
  同级，环境不再硬编码 `.conda/envs`。
- **运行底座（M2）**：`env_router.resolve_rscript` 四层优先级（模块配置 → `ULTIMATE_R_ENV_*`/
  `ULTIMATE_RSCRIPT` → `environment_root()/<env>/bin/Rscript` → `PATH`）；显式配置无效时报错
  而非静默换解释器。wheel 打包已覆盖 `scripts/R/*.R`、`config/default.yaml`、Slurm 模板与环境定义。
- **统计严谨性（M3）**：所有固定指数生存曲线已删除，改由 `survival_backend` 计算真实
  KM/Log-rank/Cox（lifelines 或 R `survival`，两者皆缺时显式 `skipped`）。Python 侧两样本检验
  改为 Welch t（`scipy.stats.ttest_ind(equal_var=False)` + Welch df + BH FDR），缺 scipy 时标注
  `normal_approx_fallback_not_for_production`。伪 WGCNA 四分位分箱已删除，改为真实
  `pickSoftThreshold` + `blockwiseModules`（样本 < 15 或缺包显式 `skipped`）。Python 不再向
  `.rds`/`.RData` 写 JSON。`report.html` 增加顶部摘要卡片、`<img>` 图表网格、Top 10 预览与
  解释边界框，并支持 `report.embed_images: true` 单文件内嵌。
- **单细胞与多组学纵深（M4）**：Harmony 批次校正（`X_pca_harmony` + `batch_integration_metrics.tsv`
  含 kNN 批次熵/iLISI 近似与 cluster 保持度 + `umap_batch_before_after.png`，缺 harmonypy 显式
  `skipped`）；配对 `edgeR` pseudobulk DE 已进入 `ultimate run` 主流水线并注册到
  `advanced_backend_execution.tsv`；CopyKAT 具备**三层闸**（显式开关 + raw integer count 门禁 +
  `thread_budget` 传 `--ncores`）；新增本地 GMT 富集 `enrichment_backend`（拒绝远程 URL，
  无 GMT 显式 `skipped`，不联网）；空间邻域统计优先真实 `squidpy.gr`（记录真实版本与
  `execution_mode=real_squidpy_gr`），缺失时降级并标注 `fallback_no_permutation_pvalue`。
- **开源上游与回归套件（M5）**：新增 `ultimate prepare-upstream --route salmon|starsolo`，生成
  `run_upstream.sbatch` / `expected_import_project.yaml` / `checksum_manifest.tsv` /
  `job_package_manifest.json` 四件产物（只生成、不自动提交），六类输入/索引问题明确 `blocked`；
  `envs/xielab.upstream.yml` 与 `slurm/xielab_setup_upstream_env.sbatch` 支持在
  `compute_pro` 上自建 fastp/salmon/STAR 工具链；`slurm/xielab_e2e_regression_suite.sbatch`
  串联 doctor → bulk → scRNA → 蛋白/生存 → validation-index/audit-production 并输出中文
  `xielab_e2e_report.md`。部署细节见 `docs/XIELAB_OPERATIONS.md`。
- **M6 集群端到端上游与全栈闭环验收（M6-P0 + M6-P1 + 全栈闭环 `PASS`）**：
  - P0：`expected_import_project.yaml` 与 `run_manifest.json` 公共 CLI 闭环及多样本独立校验、三元组 SHA-256 绑定与 4 类反例回归全部通过；
  - P1：在 `xielab compute_pro` 完成 Bulk Salmon（`213488` → `213490`，80 特征 × 4 样本）与 scRNA STARsolo（`213489` → `213491`，`SC_S1` 80×125 与 `SC_S2` 80×40）全链路实跑，并由 `ultimate check-upstream-acceptance` 输出 `status=pass`，结果记录于 `deliverables/m6/M6_ACCEPTANCE.json`；
  - **集群真环境 `pytest` 全量验证（`213751`，`node66`，耗时 `00:25:16`，退出码 `0:0`）**：`ultimate-core` + `ultimate-bio-r` **`629 passed, 25 skipped, 0 failed`**；`ultimate-scrna` + `pseudobulk-edger`（`SingleCellExperiment`/`edgeR`/`scanpy`/`anndata`/`h5py`）**`86 passed, 0 skipped, 0 failed`**。
  - **5 阶段公开基准 E2E 闭环实跑（`213750`，`node65`，耗时 `00:03:35`，退出码 `0:0`）**：在 `${ULTIMATE_RUNTIME_ROOT}/cache/public_data/` 就绪非患者公开基准集（`rnaseq_airway` 17,108 基因×8、`h.all.v7.5.1.symbols.gmt`、`kang18_8vs8.h5ad` 29,065 细胞 + `Immune_All_Low.pkl`、`UbiLength` 蛋白组 1,784×6、`survival::lung` 228 例），运行 `slurm/xielab_e2e_regression_suite.sbatch` 产出 `validations/xielab_e2e_20261128T093447Z/xielab_e2e_report.md`：**`ready=17, skipped=0, blocked=0, failed=0`**（真 R `DESeq2` + 离线 R `fgsea`、`CellTypist` + `paired_edgeR`、真 R `limma::lmFit+eBayes (R)`、真 `survival`、`report.html`、`validation-index`、`audit-production` 全通）。
  - 闭环修复（已合入主线并复测通过）：`bulk.py` 修复蛋白组 `limma` 结果列名映射（`feature_id`/`log2FC` → `protein_id`/`log2_abundance_delta`）并将局部 GMT 富集移至形式化 DE 后端之后执行；新增 `scripts/R/enrichment_fgsea_backend.R` 并更新 `enrichment_backend.py` 支持无 `gseapy` 时离线调用 R `fgsea` 与识别 `DESeq2` `log2FoldChange`/`stat`；`config.py:load_samples` 支持回落读取已启用模块的 `samplesheet`；修正单测假 `Rscript` 在有 `ULTIMATE_ENV_ROOT` 时的环境隔离。

## 下一推进点

- 上一优先级 **M1–M6 及集群真环境闭环已全部完成**；下一优先级为在 `xielab compute_pro` 上以显式开关 + 小抽样完成 `CopyKAT` 首次真实执行，再扩展到 `inferCNV`，并按需补齐 Salmon `tximport` 转录本→基因汇总链路。
- 备选方向：scATAC/Multiome fragments-level 验证、CellChat/NicheNet 通讯专项、clinical association/survival 模板。
- 真实或较重计算继续通过 Slurm 提交，不在登录节点运行。
