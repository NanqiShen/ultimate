"""Real survival analysis for the clinical_assoc / publicdb / single_gene modules.

The previous implementation drew a fixed exponential decay curve that had
nothing to do with the input data. This module computes an actual Kaplan-Meier
curve, a log-rank test and a Cox proportional-hazards ratio from the clinical
table, using ``lifelines`` when available, otherwise an R ``survival`` script
through :func:`ultimate.env_router.resolve_rscript`.

When neither engine exists, or when the clinical table carries no usable
``survival_time``/``event`` columns, the caller receives an explicit
``status=skipped`` result and a "not computed" figure; a curve is never invented.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ultimate.env_router import resolve_packaged_asset, resolve_rscript
from ultimate.plot_style import apply_clinical_journal_style, save_figure

SURVIVAL_COLUMN = "survival_time"
EVENT_COLUMN = "event"


@dataclass
class SurvivalResult:
    status: str
    reason: str = ""
    km_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    logrank_p: float | None = None
    cox_hr: float | None = None
    cox_hr_low: float | None = None
    cox_hr_high: float | None = None
    n_samples: int = 0
    engine: str = ""
    figure_path: str = ""

    def to_manifest(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "reason": self.reason,
            "logrank_p": self.logrank_p,
            "cox_hr": self.cox_hr,
            "cox_hr_95ci_low": self.cox_hr_low,
            "cox_hr_95ci_high": self.cox_hr_high,
            "n_samples": self.n_samples,
            "engine": self.engine,
            "figure": self.figure_path,
        }


def validate_survival_columns(clinical: pd.DataFrame) -> tuple[bool, str]:
    """Return ``(usable, reason)`` for the clinical survival columns."""
    if clinical is None or clinical.empty:
        return False, "missing_clinical_table"
    if SURVIVAL_COLUMN not in clinical.columns or EVENT_COLUMN not in clinical.columns:
        return False, f"missing_survival_columns:{SURVIVAL_COLUMN},{EVENT_COLUMN}"
    time_values = pd.to_numeric(clinical[SURVIVAL_COLUMN], errors="coerce")
    if time_values.isna().all() or (time_values < 0).all():
        return False, "survival_time_not_nonnegative_numeric"
    event_values = pd.to_numeric(clinical[EVENT_COLUMN], errors="coerce")
    if not set(event_values.dropna().unique()).issubset({0.0, 1.0}):
        return False, "event_column_not_binary_0_1"
    return True, ""


def _coerce_survival(clinical: pd.DataFrame) -> pd.DataFrame:
    frame = clinical.copy()
    frame[SURVIVAL_COLUMN] = pd.to_numeric(frame[SURVIVAL_COLUMN], errors="coerce")
    frame[EVENT_COLUMN] = pd.to_numeric(frame[EVENT_COLUMN], errors="coerce").fillna(0.0)
    return frame


def _feature_groups(expression: pd.Series) -> pd.Series:
    """Split samples into high/low by expression median, aligned by index."""
    threshold = expression.median()
    return pd.Series(np.where(expression >= threshold, "high", "low"), index=expression.index, name="risk_group")


def run_survival_analysis(
    clinical: pd.DataFrame,
    expression: pd.Series,
    *,
    figures_dir: Path,
    tables_dir: Path,
    label: str,
    module_name: str = "clinical_assoc",
) -> SurvivalResult:
    """Compute a real KM curve and log-rank/Cox statistics for one feature.

    ``expression`` is a per-sample series aligned with ``clinical`` on
    ``sample_id``; its index is used to join with the clinical table.
    """
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    usable, reason = validate_survival_columns(clinical)
    if not usable:
        figure = _write_skipped_figure(figures_dir, label, reason)
        return SurvivalResult(status="skipped", reason=reason, figure_path=str(figure))

    clinical_n = _coerce_survival(clinical).set_index(clinical.index)
    sample_ids = expression.index.astype(str)
    clinical_n = clinical_n.assign(_sample_id=clinical_n.get("sample_id", pd.Series(dtype=str)).astype(str))
    if "sample_id" in clinical.columns:
        clinical_n = clinical_n.set_index("sample_id")
    else:
        clinical_n = clinical_n.set_index(clinical_n["_sample_id"])

    groups = _feature_groups(expression).astype(str)
    combined = pd.DataFrame(
        {
            SURVIVAL_COLUMN: clinical_n[SURVIVAL_COLUMN],
            EVENT_COLUMN: clinical_n[EVENT_COLUMN],
            "group": groups.reindex(clinical_n.index),
        }
    ).dropna(subset=[SURVIVAL_COLUMN, "group"])
    if combined.empty or combined["group"].nunique() < 2:
        figure = _write_skipped_figure(figures_dir, label, "insufficient_samples_for_two_groups")
        return SurvivalResult(status="skipped", reason="insufficient_samples_for_two_groups", figure_path=str(figure))

    engine = _available_engine(module_name)
    if engine == "lifelines":
        return _run_lifelines(combined, label, figures_dir, tables_dir)
    if engine == "r_survival":
        return _run_r_survival(combined, label, figures_dir, tables_dir, module_name)
    figure = _write_skipped_figure(figures_dir, label, "no_survival_engine_lifelines_or_R")
    return SurvivalResult(status="skipped", reason="no_survival_engine_lifelines_or_R", figure_path=str(figure))


def _available_engine(module_name: str) -> str:
    try:
        import lifelines  # noqa: F401

        return "lifelines"
    except ImportError:
        pass
    if resolve_rscript(module_name) is not None:
        return "r_survival"
    return ""


def _run_lifelines(combined: pd.DataFrame, label: str, figures_dir: Path, tables_dir: Path) -> SurvivalResult:
    import lifelines
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test

    durations = combined[SURVIVAL_COLUMN].to_numpy(dtype=float)
    events = combined[EVENT_COLUMN].to_numpy(dtype=float)
    group = combined["group"].to_numpy()

    kmf = KaplanMeierFitter()
    km_tables: list[pd.DataFrame] = []
    for name in ("high", "low"):
        mask = group == name
        if mask.sum() == 0:
            continue
        kmf.fit(durations[mask], events[mask], label=name)
        table = kmf.survival_function_at_times(np.unique(durations[mask])).rename_axis("time").reset_index()
        table["group"] = name
        km_tables.append(table)
    km_table = pd.concat(km_tables, ignore_index=True) if km_tables else pd.DataFrame()

    high_mask = group == "high"
    low_mask = group == "low"
    logrank_p: float | None = None
    if high_mask.sum() > 0 and low_mask.sum() > 0:
        result = logrank_test(durations[high_mask], durations[low_mask], events[high_mask], events[low_mask])
        logrank_p = float(result.p_value)

    cox = combined.assign(group_high=(group == "high").astype(int))
    cox_hr = cox_hr_low = cox_hr_high = None
    if cox["group_high"].nunique() == 2 and cox[EVENT_COLUMN].nunique() >= 1:
        try:
            cph = CoxPHFitter()
            cph.fit(cox[[SURVIVAL_COLUMN, EVENT_COLUMN, "group_high"]], duration_col=SURVIVAL_COLUMN, event_col=EVENT_COLUMN)
            with np.errstate(over="ignore", invalid="ignore"):
                cox_hr = float(np.exp(cph.params_["group_high"]))
                ci = cph.confidence_intervals_.loc["group_high"]
                cox_hr_low = float(np.exp(ci.iloc[0]))
                cox_hr_high = float(np.exp(ci.iloc[1]))
            # A separated/small cohort can yield infinite bounds; report as None
            # instead of leaking inf/NaN into the delivery manifest.
            for key in ("cox_hr", "cox_hr_low", "cox_hr_high"):
                value = locals()[key]
                if not np.isfinite(value):
                    if key == "cox_hr":
                        cox_hr = None
                    elif key == "cox_hr_low":
                        cox_hr_low = None
                    else:
                        cox_hr_high = None
        except Exception:
            cox_hr = cox_hr_low = cox_hr_high = None

    figure = _plot_km(durations, events, group, label, figures_dir, logrank_p)

    if not km_table.empty:
        km_table.to_csv(tables_dir / "survival_results.tsv", sep="\t", index=False)

    return SurvivalResult(
        status="ready",
        km_table=km_table,
        logrank_p=logrank_p,
        cox_hr=cox_hr,
        cox_hr_low=cox_hr_low,
        cox_hr_high=cox_hr_high,
        n_samples=int(combined.shape[0]),
        engine="lifelines",
        figure_path=str(figure),
    )


def _run_r_survival(combined: pd.DataFrame, label: str, figures_dir: Path, tables_dir: Path, module_name: str = "clinical_assoc") -> SurvivalResult:
    """Fallback to an R survival script when lifelines is unavailable."""
    script = resolve_packaged_asset("scripts", "R", "survival_backend.R")
    input_tsv = tables_dir / "survival_input.tsv"
    output_tsv = tables_dir / "survival_results.tsv"
    combined.rename(columns={SURVIVAL_COLUMN: "time", EVENT_COLUMN: "event"}).to_csv(input_tsv, sep="\t", index=False)
    rscript = resolve_rscript(module_name)
    figure = figures_dir / "km_curve.png"
    command = [
        str(rscript),
        str(script),
        "--input",
        str(input_tsv),
        "--output",
        str(output_tsv),
        "--figure",
        str(figure),
        "--label",
        label,
    ]
    try:
        completed = subprocess.run(command, check=True, text=True, capture_output=True, timeout=300)
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        return SurvivalResult(status="skipped", reason=f"r_survival_failed:{type(exc).__name__}", figure_path=str(_write_skipped_figure(figures_dir, label, "r_survival_failed")))
    try:
        payload = json.loads(completed.stdout or "{}")
    except (json.JSONDecodeError, TypeError):
        payload = {}
    # R exits 0 even when it emits status=skipped (e.g. missing survival
    # package). Never turn that into a ready result pointing at a missing plot.
    if payload.get("status") != "ready" or not figure.is_file():
        reason = str(payload.get("reason") or "r_survival_skipped")
        skipped_figure = _write_skipped_figure(figures_dir, label, reason)
        return SurvivalResult(status="skipped", reason=reason, figure_path=str(skipped_figure))
    return SurvivalResult(
        status="ready",
        logrank_p=payload.get("logrank_p"),
        cox_hr=payload.get("cox_hr"),
        cox_hr_low=payload.get("cox_hr_low"),
        cox_hr_high=payload.get("cox_hr_high"),
        n_samples=int(combined.shape[0]),
        engine="r_survival",
        figure_path=str(figure),
    )


def _plot_km(durations: np.ndarray, events: np.ndarray, group: np.ndarray, label: str, figures_dir: Path, logrank_p: float | None) -> Path:
    from lifelines import KaplanMeierFitter

    tokens = apply_clinical_journal_style()
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6.4, 4.8))
    kmf = KaplanMeierFitter()
    for name, color_key in (("high", "case"), ("low", "control")):
        mask = group == name
        if mask.sum() == 0:
            continue
        kmf.fit(durations[mask], events[mask], label=name)
        kmf.plot_survival_function(color=tokens[color_key], linewidth=1.8, ci_show=True)
    title = f"Kaplan-Meier: {label}"
    if logrank_p is not None:
        title += f"\nLog-rank p = {logrank_p:.3g}"
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.legend(title="Expression group")
    plt.tight_layout()
    path = figures_dir / "km_curve.png"
    save_figure(path, style=tokens)
    plt.close()
    return path


def _write_skipped_figure(figures_dir: Path, label: str, reason: str) -> Path:
    tokens = apply_clinical_journal_style()
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6.4, 4.8))
    plt.text(0.5, 0.5, f"Survival analysis skipped\n{reason}", ha="center", va="center", fontsize=12)
    plt.axis("off")
    path = figures_dir / "km_curve.png"
    save_figure(path, style=tokens)
    plt.close()
    return path
