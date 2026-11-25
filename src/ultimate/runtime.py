"""Execution controls shared by CLI diagnostics and the existing pipeline."""
from __future__ import annotations

import json
import os
import sys
from contextlib import contextmanager, nullcontext
from dataclasses import asdict, dataclass
from time import monotonic
from typing import TextIO

_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
                "RCPP_PARALLEL_NUM_THREADS")


@dataclass(frozen=True)
class ThreadBudget:
    threads: int
    requested: int
    available: int
    source: str
    capped: bool

    def as_dict(self):
        return asdict(self)


def _positive(value, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be a positive integer")
    text = str(value)
    if not text.isascii() or not text.isdecimal() or int(text) < 1:
        raise ValueError(f"{field} must be a positive integer")
    return int(text)


def thread_budget(config: dict, override: int | None = None) -> ThreadBudget:
    available = os.cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        try:
            available = len(os.sched_getaffinity(0)) or available
        except OSError:
            pass
    allocation = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocation:
        available = min(available, _positive(allocation, "SLURM_CPUS_PER_TASK"))
    elif os.environ.get("SLURM_JOB_ID"):
        # Slurm defaults to one CPU per task if no explicit count was requested.
        available = min(available, 1)
    runtime = config.get("runtime", {})
    if not isinstance(runtime, dict):
        raise ValueError("runtime must be a mapping")
    if override is not None:
        value, source = override, "--threads"
    elif os.environ.get("ULTIMATE_THREADS"):
        value, source = os.environ["ULTIMATE_THREADS"], "ULTIMATE_THREADS"
    elif runtime.get("threads") is not None:
        value, source = runtime["threads"], "runtime.threads"
    else:
        value, source = (available, "slurm") if allocation else (min(available, 4), "conservative_default")
    requested = _positive(value, source)
    return ThreadBudget(min(available, requested), requested, available, source, requested > available)


@contextmanager
def apply_thread_budget(budget: ThreadBudget):
    """Limit supported native pools and inherited env; restore both on exit."""
    from threadpoolctl import threadpool_limits

    previous = {key: os.environ.get(key) for key in _THREAD_VARS}
    numexpr = sys.modules.get("numexpr")
    numexpr_previous = None
    try:
        for k in _THREAD_VARS:
            os.environ[k] = str(budget.threads)
        if numexpr is not None and hasattr(numexpr, "get_num_threads"):
            limit = min(budget.threads, numexpr.MAX_THREADS)
            numexpr_previous = numexpr.set_num_threads(limit)
        with threadpool_limits(limits=budget.threads):
            yield
    finally:
        if numexpr_previous is not None and numexpr_previous > 0:
            numexpr.set_num_threads(numexpr_previous)
        for k, v in previous.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


class RunTrace:
    """Optional JSONL; does not serialize paths, input data or exception messages."""
    def __init__(self, stream: TextIO):
        self.stream = stream

    def event(self, stage, status, **fields):
        self.stream.write(json.dumps({"stage": stage, "status": status, **fields}, sort_keys=True) + "\n")
        self.stream.flush()

    @contextmanager
    def stage(self, stage, **fields):
        start = monotonic()
        self.event(stage, "started", **fields)
        try:
            yield
        except BaseException as exc:
            self.event(stage, "failed", duration_seconds=round(monotonic()-start, 6),
                       error_type=type(exc).__name__, **fields)
            raise
        else:
            self.event(stage, "completed", duration_seconds=round(monotonic()-start, 6), **fields)


def trace_stage(trace: RunTrace | None, stage: str, **fields):
    return trace.stage(stage, **fields) if trace else nullcontext()
