"""Strict acceptance checker for an upstream -> import -> run chain.

The M5 five-stage suite treated ``blocked``/``skipped`` stages as acceptable and
only failed on ``failed``, which meant "nothing actually ran" could still exit 0.
This checker inverts that: every *required* check must be present **and** pass.
A missing, skipped, blocked or failed required check makes the whole run fail.

It is intentionally data-driven: callers build a :class:`Check` list, and the
checker only decides pass/fail plus writes a machine-readable report.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

ACCEPTANCE_SCHEMA = "ultimate.m6.acceptance.v1"

STATUS_PASS = "pass"
STATUS_FAIL = "fail"
STATUS_SKIPPED = "skipped"
STATUS_BLOCKED = "blocked"
STATUS_MISSING = "missing"

# Any status other than `pass` fails a required check.
PASSING_STATUSES = frozenset({STATUS_PASS})


@dataclass
class Check:
    name: str
    status: str
    detail: str = ""
    required: bool = True

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "detail": self.detail,
            "required": self.required,
        }


@dataclass
class RouteAcceptance:
    route: str
    checks: list[Check] = field(default_factory=list)

    def add(self, name: str, status: str, detail: str = "", required: bool = True) -> Check:
        check = Check(name=name, status=status, detail=detail, required=required)
        self.checks.append(check)
        return check

    @property
    def failures(self) -> list[Check]:
        return [check for check in self.checks if check.required and check.status not in PASSING_STATUSES]

    def to_dict(self) -> dict[str, Any]:
        return {
            "route": self.route,
            "status": STATUS_PASS if not self.failures else STATUS_FAIL,
            "checks": [check.to_dict() for check in self.checks],
            "failed_checks": [check.name for check in self.failures],
        }


def check_expected_files(acceptance: RouteAcceptance, prefix: str, paths: Iterable[Path]) -> None:
    """Every expected artefact must exist and be non-empty."""
    for path in paths:
        path = Path(path)
        if not path.exists():
            acceptance.add(f"{prefix}:exists:{path.name}", STATUS_MISSING, str(path))
        elif path.is_file() and path.stat().st_size == 0:
            acceptance.add(f"{prefix}:nonempty:{path.name}", STATUS_FAIL, f"empty file: {path}")
        else:
            acceptance.add(f"{prefix}:exists:{path.name}", STATUS_PASS, str(path))


def check_matrix_contract(
    acceptance: RouteAcceptance,
    prefix: str,
    matrix_path: Path,
    *,
    expected_samples: list[str] | None = None,
    require_integer: bool = False,
    require_nonnegative: bool = True,
    require_nonzero_signal: bool = True,
) -> dict[str, Any]:
    """Validate a counts matrix's shape, samples, uniqueness and value domain."""
    import pandas as pd

    matrix_path = Path(matrix_path)
    if not matrix_path.is_file():
        acceptance.add(f"{prefix}:matrix_readable", STATUS_MISSING, str(matrix_path))
        return {}
    try:
        frame = pd.read_csv(matrix_path, sep="\t")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        acceptance.add(f"{prefix}:matrix_parse", STATUS_FAIL, f"{type(exc).__name__}: {exc}")
        return {}
    if frame.shape[0] == 0 or frame.shape[1] < 2:
        acceptance.add(f"{prefix}:matrix_shape", STATUS_FAIL, f"shape={frame.shape}")
        return {"shape": list(frame.shape)}

    index_name = frame.columns[0]
    feature_ids = frame[index_name].astype(str)
    acceptance.add(
        f"{prefix}:feature_ids_unique",
        STATUS_PASS if feature_ids.is_unique else STATUS_FAIL,
        f"{feature_ids.nunique()}/{len(feature_ids)} unique",
    )
    sample_columns = [str(col) for col in frame.columns[1:]]
    if expected_samples is not None:
        missing = sorted(set(expected_samples) - set(sample_columns))
        extra = sorted(set(sample_columns) - set(expected_samples))
        acceptance.add(
            f"{prefix}:samples_match_manifest",
            STATUS_PASS if not missing and not extra else STATUS_FAIL,
            f"missing={missing} extra={extra}",
        )
    values = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    acceptance.add(
        f"{prefix}:values_finite",
        STATUS_PASS if not values.isna().to_numpy().any() else STATUS_FAIL,
        f"nan_cells={int(values.isna().to_numpy().sum())}",
    )
    if require_nonnegative:
        negative = int((values.fillna(0) < 0).to_numpy().sum())
        acceptance.add(f"{prefix}:values_nonnegative", STATUS_PASS if negative == 0 else STATUS_FAIL, f"negative={negative}")
    if require_integer:
        import numpy as np

        non_integer = int((~np.isclose(values.fillna(0).to_numpy(), np.round(values.fillna(0).to_numpy()))).sum())
        acceptance.add(f"{prefix}:values_integer", STATUS_PASS if non_integer == 0 else STATUS_FAIL, f"non_integer={non_integer}")
    if require_nonzero_signal:
        total = float(values.fillna(0).to_numpy().sum())
        acceptance.add(f"{prefix}:nonzero_signal", STATUS_PASS if total > 0 else STATUS_FAIL, f"total={total}")
    return {
        "shape": list(frame.shape),
        "features": int(len(feature_ids)),
        "samples": sample_columns,
        "index_column": index_name,
    }


def check_10x_directory(acceptance: RouteAcceptance, prefix: str, directory: Path) -> dict[str, Any]:
    """A 10x ``filtered/`` directory must contain matrix/barcodes/features."""
    import gzip
    import io

    directory = Path(directory)
    info: dict[str, Any] = {}
    for stem in ("matrix", "barcodes", "features"):
        candidates = [directory / f"{stem}.mtx", directory / f"{stem}.mtx.gz",
                      directory / f"{stem}.tsv", directory / f"{stem}.tsv.gz"]
        found = next((c for c in candidates if c.is_file()), None)
        acceptance.add(
            f"{prefix}:10x:{stem}",
            STATUS_PASS if found else STATUS_MISSING,
            str(found or directory / stem),
        )
        if found is not None:
            info[stem] = found
    # Cross-check declared dimensions against the actual barcode/feature counts.
    if {"matrix", "barcodes", "features"} <= set(info):
        def _lines(path: Path) -> list[str]:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
                return [line for line in handle.read().splitlines() if line.strip()]
        barcodes = _lines(info["barcodes"])
        features = _lines(info["features"])
        matrix_lines = [line for line in _lines(info["matrix"]) if not line.startswith("%")]
        try:
            n_features, n_barcodes, _nnz = (int(x) for x in matrix_lines[0].split()[:3])
        except (ValueError, IndexError):
            acceptance.add(f"{prefix}:10x:shape_header", STATUS_FAIL, "unparsable matrix header")
        else:
            acceptance.add(
                f"{prefix}:10x:dims_consistent",
                STATUS_PASS if n_features == len(features) and n_barcodes == len(barcodes) else STATUS_FAIL,
                f"matrix={n_features}x{n_barcodes} features={len(features)} barcodes={len(barcodes)}",
            )
            info["shape"] = [n_features, n_barcodes]
    return info


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_downstream_consumed(
    acceptance: RouteAcceptance,
    prefix: str,
    run_manifest_path: Path,
    *,
    expected_input: Path,
    min_cells: int = 1,
) -> dict[str, Any]:
    """Prove Ultimate actually read *this* upstream artefact in this run.

    Reads the downstream run manifest and requires that its recorded input path
    resolves to the artefact this suite just produced (not a demo fallback and
    not a pre-existing matrix).
    """
    run_manifest_path = Path(run_manifest_path)
    if not run_manifest_path.is_file():
        acceptance.add(f"{prefix}:downstream_manifest", STATUS_MISSING, str(run_manifest_path))
        return {}
    try:
        payload = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        acceptance.add(f"{prefix}:downstream_manifest", STATUS_FAIL, f"{type(exc).__name__}")
        return {}
    recorded = str(payload.get("input_path") or payload.get("input_matrix") or "")
    resolved = Path(recorded).resolve() if recorded else None
    expected = Path(expected_input).resolve()
    acceptance.add(
        f"{prefix}:downstream_consumed_this_artefact",
        STATUS_PASS if resolved == expected else STATUS_FAIL,
        f"manifest={recorded} expected={expected}",
    )
    status = str(payload.get("status") or "")
    acceptance.add(
        f"{prefix}:downstream_status",
        STATUS_PASS if status in {"ready", "complete_python_bulk_backend"} else STATUS_FAIL,
        f"status={status}",
    )
    delivery_allowed = payload.get("delivery_allowed")
    acceptance.add(
        f"{prefix}:delivery_allowed_false",
        STATUS_PASS if delivery_allowed is not True else STATUS_FAIL,
        f"delivery_allowed={delivery_allowed} (false is the expected governance result)",
        required=True,
    )
    n_cells = payload.get("n_cells") or payload.get("n_samples")
    if isinstance(n_cells, int):
        acceptance.add(
            f"{prefix}:downstream_nonempty",
            STATUS_PASS if n_cells >= min_cells else STATUS_FAIL,
            f"n={n_cells}",
        )
    return {"recorded_input": recorded, "status": status, "n": n_cells}


def check_sample_count_agreement(
    acceptance: RouteAcceptance,
    prefix: str,
    matrix_path: Path,
    samplesheet_path: Path,
) -> list[str]:
    """Matrix sample columns must equal the normalised samplesheet sample ids."""
    import pandas as pd

    try:
        sheet = pd.read_csv(samplesheet_path, sep="\t")
        frame = pd.read_csv(matrix_path, sep="\t", nrows=0)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        acceptance.add(f"{prefix}:sample_count_agreement", STATUS_FAIL, "unreadable inputs")
        return []
    expected = [str(x) for x in sheet["sample_id"].tolist()]
    actual = [str(col) for col in frame.columns[1:]]
    acceptance.add(
        f"{prefix}:sample_count_agreement",
        STATUS_PASS if expected == actual else STATUS_FAIL,
        f"sheet={expected} matrix={actual}",
    )
    return expected


def evaluate(routes: list[RouteAcceptance]) -> dict[str, Any]:
    """Aggregate route results; overall status is pass only when all pass."""
    failed = [route.route for route in routes if route.failures]
    return {
        "schema": ACCEPTANCE_SCHEMA,
        "status": STATUS_PASS if not failed else STATUS_FAIL,
        "failed_routes": failed,
        "routes": [route.to_dict() for route in routes],
    }


def write_acceptance(path: Path, payload: dict[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def exit_code(payload: dict[str, Any]) -> int:
    """0 only when the aggregate status is pass; anything else is non-zero."""
    return 0 if payload.get("status") == STATUS_PASS else 1
