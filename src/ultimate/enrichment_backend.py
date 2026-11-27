"""Local-only pathway enrichment (GSEA prerank / over-representation).

XieLab runs offline: this module never contacts Enrichr, KEGG or any remote
service. Gene sets must come from a local GMT file (``resources.msigdb``,
``modules.<module>.gmt_path`` or ``@resources/msigdb/*.gmt``). When no GMT file
or no enrichment engine is available, the backend records an explicit
``status=skipped`` with the concrete reason instead of downloading anything.
"""

from __future__ import annotations

import glob
import importlib.util
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ultimate.plot_style import apply_clinical_journal_style, save_figure

ENRICHMENT_BACKEND_ID = "enrichment.gsea_ora_local"
ENRICHMENT_WARNING = (
    "通路富集结果是基于本地基因集文件的统计富集，属于假设生成证据；"
    "必须结合基因集版本、背景基因数与被检验基因数量解释，不能写成机制证明。"
)
MIN_RANKED_GENES = 10


def _module_available(name: str) -> bool:
    """True when a module can actually be imported.

    ``importlib.util.find_spec`` raises ``ValueError`` for partially initialised
    modules (``__spec__ is None``); a probe must never crash the backend.
    """
    import importlib
    import importlib.util

    try:
        if importlib.util.find_spec(name) is not None:
            return True
    except (ImportError, ValueError, AttributeError):
        pass
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False




def _resolve_gmt_paths(config: dict[str, Any], module_cfg: dict[str, Any], module_name: str) -> list[Path]:
    """Collect local GMT files; never a URL."""
    candidates: list[str] = []
    for value in (
        module_cfg.get("gmt_path"),
        module_cfg.get("gmt"),
        (module_cfg.get("enrichment") or {}).get("gmt_path") if isinstance(module_cfg.get("enrichment"), dict) else None,
    ):
        if value:
            candidates.append(str(value))
    resources = config.get("resources") if isinstance(config.get("resources"), dict) else {}
    if resources.get("msigdb"):
        candidates.append(str(resources["msigdb"]))
    paths: list[Path] = []
    for candidate in candidates:
        if candidate.startswith(("http://", "https://", "ftp://")):
            # Explicitly refuse remote gene sets: offline cluster policy.
            continue
        text = candidate.replace("@resources/", "").replace("@project/", "")
        path = Path(text).expanduser()
        if path.is_file() and path.suffix.lower() == ".gmt":
            paths.append(path)
        elif path.is_dir():
            paths.extend(sorted(Path(p) for p in glob.glob(str(path / "*.gmt"))))
        else:
            sibling = path.with_suffix(".gmt")
            if sibling.is_file():
                paths.append(sibling)
    unique: list[Path] = []
    for path in paths:
        if path not in unique:
            unique.append(path)
    return unique


def _read_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        genes = [gene.strip() for gene in parts[2:] if gene.strip()]
        if genes:
            gene_sets[parts[0].strip()] = genes
    return gene_sets


def _ranked_genes(de_table: pd.DataFrame, *, min_genes: int = MIN_RANKED_GENES) -> pd.Series:
    """Ranked gene -> score series from a DE table (log2FC ordered, padj available)."""
    if de_table is None or de_table.empty:
        return pd.Series(dtype=float)
    feature_column = next(
        (col for col in ("feature_id", "gene", "gene_id", "protein_id", "region_id", "gene_symbol") if col in de_table.columns),
        None,
    )
    if feature_column is None:
        return pd.Series(dtype=float)
    score_column = next(
        (col for col in ("log2FC", "log2_abundance_delta", "m_value_delta", "statistic", "score", "logfoldchanges") if col in de_table.columns),
        None,
    )
    if score_column is None:
        return pd.Series(dtype=float)
    frame = de_table[[feature_column, score_column]].copy()
    frame[score_column] = pd.to_numeric(frame[score_column], errors="coerce")
    frame = frame.dropna().drop_duplicates(subset=[feature_column])
    if len(frame) < min_genes:
        return pd.Series(dtype=float)
    return pd.Series(frame[score_column].to_numpy(), index=frame[feature_column].astype(str).to_numpy()).sort_values(ascending=False)


def _write_status(path: Path, status: str, reason: str, analysis_level: str, **extra: Any) -> None:
    row = {
        "backend_id": ENRICHMENT_BACKEND_ID,
        "status": status,
        "reason": reason,
        "analysis_level": analysis_level,
        "method_boundary": ENRICHMENT_WARNING,
        "delivery_allowed": False,
        "network_access": "none_local_gmt_only",
    }
    row.update(extra)
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)


def _skip(
    *,
    tables_dir: Path,
    status_path: Path,
    reason: str,
    analysis_level: str,
    artifacts: dict[str, str],
) -> dict[str, Any]:
    _write_status(status_path, "skipped", reason, analysis_level)
    result_path = tables_dir / "enrichment_results.tsv"
    pd.DataFrame(
        [{"backend_id": ENRICHMENT_BACKEND_ID, "status": "skipped", "reason": reason, "analysis_level": analysis_level, "delivery_allowed": False}]
    ).to_csv(result_path, sep="\t", index=False)
    return {"status": "skipped", "reason": reason, "artifacts": artifacts}


def run_enrichment(
    *,
    de_table: pd.DataFrame,
    config: dict[str, Any],
    module_name: str,
    output_dir: Path,
    analysis_level: str = "smoke_backend",
    threads: int | None = None,
) -> dict[str, Any]:
    """Run local GSEA prerank (preferred) or ORA over a DE table."""
    tables_dir = output_dir / "results" / "tables" / module_name
    figures_dir = output_dir / "results" / "figures" / module_name
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    status_path = tables_dir / "enrichment_backend_status.tsv"
    dotplot_path = figures_dir / "enrichment_dotplot.png"
    curve_path = figures_dir / "enrichment_gsea_curve.png"
    artifacts = {
        "enrichment_backend_status": str(status_path),
        "enrichment_dotplot": str(dotplot_path),
        "enrichment_gsea_curve": str(curve_path),
    }
    modules = config.get("modules") if isinstance(config.get("modules"), dict) else {}
    module_cfg = modules.get(module_name) if isinstance(modules.get(module_name), dict) else {}

    gmt_paths = _resolve_gmt_paths(config, module_cfg or {}, module_name)
    if not gmt_paths:
        return _skip(
            tables_dir=tables_dir, status_path=status_path,
            reason="missing_local_gmt:no_network_download_allowed",
            analysis_level=analysis_level, artifacts=artifacts,
        )
    gene_sets: dict[str, list[str]] = {}
    for gmt in gmt_paths:
        try:
            gene_sets.update(_read_gmt(gmt))
        except OSError:
            continue
    if not gene_sets:
        return _skip(
            tables_dir=tables_dir, status_path=status_path,
            reason="empty_gmt_gene_sets", analysis_level=analysis_level, artifacts=artifacts,
        )

    ranked = _ranked_genes(de_table)
    if ranked.empty:
        return _skip(
            tables_dir=tables_dir, status_path=status_path,
            reason=f"insufficient_ranked_genes:<{MIN_RANKED_GENES}",
            analysis_level=analysis_level, artifacts=artifacts,
        )
    if not _module_available("gseapy"):
        return _skip(
            tables_dir=tables_dir, status_path=status_path,
            reason="dependency_missing:gseapy", analysis_level=analysis_level, artifacts=artifacts,
        )

    import gseapy as gp

    try:
        with _thread_budget(config, threads):
            prerank = gp.prerank(
                rnk=ranked,
                gene_sets=gene_sets,
                outdir=None,
                permutation_num=100,
                min_size=5,
                max_size=500,
                seed=123,
                threads=1,
                verbose=False,
                no_plot=True,
            )
    except Exception as exc:
        return _skip(
            tables_dir=tables_dir, status_path=status_path,
            reason=f"enrichment_failed:{type(exc).__name__}", analysis_level=analysis_level, artifacts=artifacts,
        )

    result = prerank.res2d.copy()
    result.insert(0, "backend_id", ENRICHMENT_BACKEND_ID)
    result["analysis_level"] = analysis_level
    result["gene_set_source"] = ";".join(str(path) for path in gmt_paths)
    result_path = tables_dir / "enrichment_results.tsv"
    result.to_csv(result_path, sep="\t", index=False)
    artifacts["enrichment_results"] = str(result_path)
    _write_status(status_path, "ready", "", analysis_level, n_gene_sets=len(gene_sets), n_ranked_genes=len(ranked))
    _plot_dotplot(result, dotplot_path)
    _plot_gsea_curve(prerank, curve_path)
    return {
        "status": "ready",
        "reason": "",
        "artifacts": artifacts,
        "n_gene_sets": len(gene_sets),
        "n_ranked_genes": len(ranked),
    }


class _thread_budget:
    """Context wrapper so enrichment honours the run-wide CPU budget."""

    def __init__(self, config: dict[str, Any], threads: int | None):
        self._config = config
        self._threads = threads
        self._ctx = None

    def __enter__(self):
        from ultimate.runtime import apply_thread_budget, thread_budget

        self._ctx = apply_thread_budget(thread_budget(self._config, self._threads))
        return self._ctx.__enter__()

    def __exit__(self, exc_type, exc, tb):
        return self._ctx.__exit__(exc_type, exc, tb)


def _plot_dotplot(result: pd.DataFrame, path: Path) -> None:
    tokens = apply_clinical_journal_style()
    frame = result.head(15).copy()
    if frame.empty or "Term" not in frame.columns:
        return
    score_col = next((c for c in ("NES", "nes", "combined_score") if c in frame.columns), None)
    if score_col is None:
        return
    scores = pd.to_numeric(frame[score_col], errors="coerce").fillna(0.0)
    fdr_col = next((c for c in ("FDR q-val", "fdr", "qvalue", "Adjusted P-value") if c in frame.columns), None)
    fdr = pd.to_numeric(frame[fdr_col], errors="coerce").fillna(1.0) if fdr_col else pd.Series(np.ones(len(frame)))
    labels = frame["Term"].astype(str).str.slice(0, 60)
    fig, ax = plt.subplots(figsize=(7.6, max(3.2, 0.32 * len(frame) + 1.4)))
    ax.scatter(scores, range(len(frame)), s=80, c=-np.log10(np.clip(fdr, 1e-10, 1)), cmap="viridis")
    ax.set_yticks(range(len(frame)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0.0, color=tokens["axis"], linewidth=0.8)
    ax.set_xlabel(str(score_col))
    ax.set_title("Local GMT enrichment (top 15)", fontsize=11)
    fig.tight_layout()
    save_figure(path, close=True)
    plt.close(fig)


def _plot_gsea_curve(prerank: Any, path: Path) -> None:
    tokens = apply_clinical_journal_style()
    ranking = getattr(prerank, "ranking", None)
    results = getattr(prerank, "results", None)
    if not isinstance(results, dict) or not results or ranking is None:
        return
    first = next(iter(results.values()))
    term_name = str(first.get("Term") or "top path")
    es = np.asarray(first.get("RES", []), dtype=float) if first.get("RES") is not None else np.array([])
    if es.size == 0:
        return
    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    ax.plot(range(es.size), es, color=tokens["case"], linewidth=1.4)
    ax.axhline(0.0, color=tokens["axis"], linewidth=0.7)
    ax.fill_between(range(es.size), 0, es, where=es > 0, color=tokens["case"], alpha=0.20)
    ax.set_title(f"GSEA running enrichment score\n{term_name} (local GMT)", fontsize=10)
    ax.set_xlabel("Ranked gene position")
    ax.set_ylabel("ES")
    fig.tight_layout()
    save_figure(path, close=True)
    plt.close(fig)
