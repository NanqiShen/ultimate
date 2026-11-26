"""Strict acceptance checker for an upstream -> import -> run chain.

The M5 five-stage suite treated ``blocked``/``skipped`` stages as acceptable and
only failed on ``failed``, which meant "nothing actually ran" could still exit 0.
This checker inverts that: every *required* check must be present **and** pass.
A missing, skipped, blocked or failed required check makes the whole run fail.

It is intentionally data-driven: callers build a :class:`Check` list, and the
checker only decides pass/fail plus writes a machine-readable report.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

ACCEPTANCE_SCHEMA = "ultimate.m6.acceptance.v1"

STATUS_PASS = "pass"
STATUS_FAIL = "fail"
STATUS_SKIPPED = "skipped"
STATUS_BLOCKED = "blocked"
STATUS_MISSING = "missing"

# Any status other than `pass` fails a required check.
PASSING_STATUSES = frozenset({STATUS_PASS})
REQUIRED_M6_ROUTES = ("salmon", "starsolo")


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
    import numpy as np
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
    raw_vals = values.to_numpy(dtype=float)
    is_finite = bool(np.isfinite(raw_vals).all())
    acceptance.add(
        f"{prefix}:values_finite",
        STATUS_PASS if is_finite else STATUS_FAIL,
        f"finite={is_finite}",
    )
    if require_nonnegative:
        negative = int((values.fillna(0) < 0).to_numpy().sum())
        acceptance.add(f"{prefix}:values_nonnegative", STATUS_PASS if negative == 0 else STATUS_FAIL, f"negative={negative}")
    if require_integer:
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


def check_salmon_quant_reconciliation(
    acceptance: RouteAcceptance,
    prefix: str,
    matrix_path: Path,
    artifacts_dir: Path,
    expected_samples: list[str],
) -> None:
    """Every sample must have a valid quant.sf whose rounded NumReads matches counts_raw.tsv."""
    import pandas as pd

    matrix_path = Path(matrix_path)
    if not matrix_path.is_file():
        acceptance.add(f"{prefix}:quant_reconciliation:matrix_missing", STATUS_MISSING, str(matrix_path))
        return

    try:
        counts_df = pd.read_csv(matrix_path, sep="\t")
    except Exception as exc:
        acceptance.add(f"{prefix}:quant_reconciliation:matrix_parse_fail", STATUS_FAIL, str(exc))
        return

    feature_col = counts_df.columns[0]
    counts_df[feature_col] = counts_df[feature_col].astype(str)
    counts_map = counts_df.set_index(feature_col)

    quant_base = artifacts_dir / "quant"
    for sample in expected_samples:
        sf_path = quant_base / sample / "quant.sf"
        if not sf_path.is_file() or sf_path.stat().st_size == 0:
            acceptance.add(
                f"{prefix}:quant_sf_present:{sample}",
                STATUS_FAIL,
                f"missing or empty quant.sf at {sf_path}",
            )
            continue
        acceptance.add(f"{prefix}:quant_sf_present:{sample}", STATUS_PASS, str(sf_path))

        if sample not in counts_map.columns:
            acceptance.add(
                f"{prefix}:quant_reconciliation:{sample}",
                STATUS_FAIL,
                f"sample {sample} not found in matrix columns",
            )
            continue

        try:
            sf_df = pd.read_csv(sf_path, sep="\t")
        except Exception as exc:
            acceptance.add(f"{prefix}:quant_sf_parse:{sample}", STATUS_FAIL, str(exc))
            continue

        if "Name" not in sf_df.columns or "NumReads" not in sf_df.columns:
            acceptance.add(
                f"{prefix}:quant_sf_columns:{sample}",
                STATUS_FAIL,
                f"missing Name or NumReads in {sf_path}",
            )
            continue

        num_reads = pd.to_numeric(sf_df["NumReads"], errors="coerce")
        if num_reads.isna().any() or not np.isfinite(num_reads.to_numpy(dtype=float)).all():
            acceptance.add(
                f"{prefix}:quant_values_finite:{sample}",
                STATUS_FAIL,
                f"non-finite values in {sf_path}",
            )
            continue

        sf_df["Name"] = sf_df["Name"].astype(str)
        sf_lookup = dict(zip(sf_df["Name"], np.round(num_reads.to_numpy(dtype=float)).astype(int)))

        # Compare rounded NumReads to counts_raw.tsv
        matrix_vals = counts_map[sample]
        mismatches = 0
        missing_features = 0
        nonfinite_in_matrix = 0
        for feat, val in matrix_vals.items():
            feat_str = str(feat)
            if feat_str not in sf_lookup:
                missing_features += 1
                continue
            try:
                fval = float(val)
                if not math.isfinite(fval):
                    nonfinite_in_matrix += 1
                    continue
                if int(round(fval)) != int(sf_lookup[feat_str]):
                    mismatches += 1
            except (ValueError, TypeError):
                nonfinite_in_matrix += 1

        if mismatches > 0 or missing_features > 0 or nonfinite_in_matrix > 0:
            acceptance.add(
                f"{prefix}:quant_reconciliation:{sample}",
                STATUS_FAIL,
                f"mismatches={mismatches} missing_features={missing_features} nonfinite={nonfinite_in_matrix} in {sample}",
            )
        else:
            acceptance.add(
                f"{prefix}:quant_reconciliation:{sample}",
                STATUS_PASS,
                f"{len(matrix_vals)} features verified with {sf_path.name}",
            )


def check_10x_directory(acceptance: RouteAcceptance, prefix: str, directory: Path) -> dict[str, Any]:
    """A 10x ``filtered/`` directory must contain matrix/barcodes/features with valid entries."""
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

    if {"matrix", "barcodes", "features"} <= set(info):
        def _lines(path: Path) -> list[str]:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
                return [line for line in handle.read().splitlines() if line.strip()]

        barcodes = _lines(info["barcodes"])
        features = _lines(info["features"])
        matrix_lines = [line for line in _lines(info["matrix"]) if not line.startswith("%")]

        if not matrix_lines:
            acceptance.add(f"{prefix}:10x:shape_header", STATUS_FAIL, "empty matrix market file")
            return info

        header = matrix_lines[0].split()
        try:
            n_features, n_barcodes, nnz = int(header[0]), int(header[1]), int(header[2])
        except (ValueError, IndexError):
            acceptance.add(f"{prefix}:10x:shape_header", STATUS_FAIL, f"unparsable header: {matrix_lines[0]}")
            return info

        acceptance.add(
            f"{prefix}:10x:dims_consistent",
            STATUS_PASS if n_features == len(features) and n_barcodes == len(barcodes) else STATUS_FAIL,
            f"matrix={n_features}x{n_barcodes} features={len(features)} barcodes={len(barcodes)}",
        )
        info["shape"] = [n_features, n_barcodes]

        # Inspect actual matrix entries (not just the header)
        data_lines = matrix_lines[1:]
        acceptance.add(
            f"{prefix}:10x:entry_count_matches_nnz",
            STATUS_PASS if len(data_lines) == nnz else STATUS_FAIL,
            f"declared_nnz={nnz} found_lines={len(data_lines)}",
        )

        entries_valid = True
        has_positive = False
        invalid_reason = ""
        for idx, line in enumerate(data_lines, start=2):
            parts = line.split()
            if len(parts) < 3:
                entries_valid = False
                invalid_reason = f"line {idx} fewer than 3 columns: {line}"
                break
            try:
                r, c = int(parts[0]), int(parts[1])
                val = float(parts[2])
            except ValueError:
                entries_valid = False
                invalid_reason = f"line {idx} non-numeric: {line}"
                break

            if r < 1 or r > n_features or c < 1 or c > n_barcodes:
                entries_valid = False
                invalid_reason = f"line {idx} coordinates out of bounds: row={r} col={c} (max={n_features},{n_barcodes})"
                break
            if not math.isfinite(val) or val < 0 or not val.is_integer():
                entries_valid = False
                invalid_reason = f"line {idx} count invalid (must be finite non-negative integer): {val}"
                break
            if val > 0:
                has_positive = True

        acceptance.add(
            f"{prefix}:10x:values_integer_nonnegative",
            STATUS_PASS if entries_valid else STATUS_FAIL,
            invalid_reason or f"{len(data_lines)} entries verified finite non-negative integer",
        )
        acceptance.add(
            f"{prefix}:10x:nonzero_signal",
            STATUS_PASS if (entries_valid and has_positive) else STATUS_FAIL,
            "positive counts found" if has_positive else "all matrix counts are zero",
        )
    return info


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _collect_manifest_checksums(payload: dict[str, Any]) -> dict[str, str]:
    """Collect any recorded input checksums from the downstream manifest."""
    checksums: dict[str, str] = {}
    for key in ("input_checksum", "input_sha256", "checksum"):
        val = payload.get(key)
        if isinstance(val, str) and len(val) == 64:
            checksums["primary"] = val

    raw_checks = payload.get("input_checksums") or payload.get("checksums")
    if isinstance(raw_checks, dict):
        for k, v in raw_checks.items():
            if isinstance(v, str) and len(v) == 64:
                checksums[Path(str(k)).name] = v
                checksums[str(k)] = v
    elif isinstance(raw_checks, list):
        for item in raw_checks:
            if isinstance(item, dict):
                p = item.get("path") or item.get("file")
                s = item.get("sha256") or item.get("checksum")
                if p and s and len(str(s)) == 64:
                    checksums[Path(str(p)).name] = str(s)
                    checksums[str(p)] = str(s)

    raw_qc = payload.get("raw_qc")
    if isinstance(raw_qc, dict):
        c = raw_qc.get("input_checksum") or raw_qc.get("checksum")
        if isinstance(c, str) and len(c) == 64:
            checksums["raw_qc"] = c
    return checksums


def check_downstream_consumed(
    acceptance: RouteAcceptance,
    prefix: str,
    run_manifest_path: Path,
    *,
    expected_inputs: list[Path] | Path,
    expected_sample_count: int | None = None,
) -> dict[str, Any]:
    """Prove Ultimate actually read *this* upstream artefact in this run with runtime evidence."""
    run_manifest_path = Path(run_manifest_path)
    if not run_manifest_path.is_file():
        acceptance.add(f"{prefix}:downstream_manifest", STATUS_MISSING, str(run_manifest_path))
        return {}
    try:
        payload = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        acceptance.add(f"{prefix}:downstream_manifest", STATUS_FAIL, f"{type(exc).__name__}")
        return {}

    # 1. Runtime identity binding
    run_id = payload.get("run_id")
    generated_at = payload.get("generated_at")
    has_runtime_id = bool(run_id and isinstance(run_id, str) and generated_at and isinstance(generated_at, str))
    acceptance.add(
        f"{prefix}:downstream_runtime_binding",
        STATUS_PASS if has_runtime_id else STATUS_FAIL,
        f"run_id={run_id} generated_at={generated_at}",
    )

    # 2. Check input path consumption
    if not isinstance(expected_inputs, list):
        expected_inputs = [expected_inputs]

    recorded_paths: list[str] = []
    for k in ("input_path", "input_matrix"):
        if payload.get(k):
            recorded_paths.append(str(payload[k]))
    for mod in payload.get("modules") or []:
        if isinstance(mod, dict):
            for k in ("input_path", "input_matrix", "matrix_path"):
                if mod.get(k):
                    recorded_paths.append(str(mod[k]))
    for samp in payload.get("samples") or []:
        if isinstance(samp, dict) and samp.get("input_path"):
            recorded_paths.append(str(samp["input_path"]))

    resolved_recorded = {str(Path(p).resolve()) for p in recorded_paths if p}
    all_consumed = True
    missing_consumptions = []
    for exp in expected_inputs:
        exp_resolved = str(Path(exp).resolve())
        if exp_resolved not in resolved_recorded:
            all_consumed = False
            missing_consumptions.append(f"{exp.parent.name}/{exp.name}" if exp.parent else exp.name)

    acceptance.add(
        f"{prefix}:downstream_consumed_this_artefact",
        STATUS_PASS if all_consumed else STATUS_FAIL,
        f"missing_consumptions={missing_consumptions} recorded={sorted(resolved_recorded)}",
    )

    # 3. Input checksum verification for all expected inputs
    manifest_checksums = _collect_manifest_checksums(payload)
    missing_checksums = []
    for exp in expected_inputs:
        target_file = exp
        if exp.is_dir():
            found_stem = None
            for stem in ("matrix.mtx.gz", "matrix.mtx"):
                if (exp / stem).is_file():
                    target_file = exp / stem
                    found_stem = stem
                    break
            if not found_stem:
                missing_checksums.append(f"missing_matrix_in:{exp.name}")
                continue
        if target_file.is_file():
            actual_sha = sha256_file(target_file)
            matched = (
                actual_sha in manifest_checksums.values()
                or manifest_checksums.get(target_file.name) == actual_sha
                or manifest_checksums.get(str(target_file)) == actual_sha
                or manifest_checksums.get("primary") == actual_sha
            )
            if not matched:
                missing_checksums.append(target_file.name)
        else:
            missing_checksums.append(f"file_not_found:{target_file.name}")

    has_all_checksums = bool(expected_inputs and not missing_checksums)
    acceptance.add(
        f"{prefix}:downstream_input_checksum",
        STATUS_PASS if has_all_checksums else STATUS_FAIL,
        f"verified {len(expected_inputs)} input checksums" if has_all_checksums else f"unmatched checksums: {missing_checksums}",
    )

    # 4. Actual dimensions
    n_cells = payload.get("n_cells")
    n_samples = payload.get("n_samples")
    n_genes = payload.get("n_genes") or payload.get("n_features")
    shape = payload.get("matrix_shape") or payload.get("dimensions")

    has_dims = False
    if isinstance(n_cells, int) and n_cells > 0 and isinstance(n_genes, int) and n_genes > 0:
        has_dims = True
    elif isinstance(n_samples, int) and n_samples > 0 and isinstance(n_genes, int) and n_genes > 0:
        has_dims = True
    elif isinstance(n_cells, int) and n_cells > 0:
        has_dims = True
    elif isinstance(shape, (list, tuple)) and len(shape) == 2 and all(isinstance(x, int) and x > 0 for x in shape):
        has_dims = True

    acceptance.add(
        f"{prefix}:downstream_actual_dimensions",
        STATUS_PASS if has_dims else STATUS_FAIL,
        f"n_cells={n_cells} n_samples={n_samples} n_genes={n_genes} shape={shape}",
    )

    # 5. Delivery allowed must be strictly boolean False
    delivery_allowed = payload.get("delivery_allowed")
    is_boolean_false = (isinstance(delivery_allowed, bool) and delivery_allowed is False)
    acceptance.add(
        f"{prefix}:delivery_allowed_boolean_false",
        STATUS_PASS if is_boolean_false else STATUS_FAIL,
        f"delivery_allowed={delivery_allowed!r} (must be boolean False, not string or omitted)",
    )

    status = str(payload.get("status") or "")
    acceptance.add(
        f"{prefix}:downstream_status",
        STATUS_PASS if status in {"ready", "complete_python_bulk_backend"} else STATUS_FAIL,
        f"status={status}",
    )
    return {"run_id": run_id, "status": status, "n_cells": n_cells, "n_samples": n_samples}


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


def evaluate(
    routes: list[RouteAcceptance],
    *,
    required_routes: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Aggregate route results; overall status is pass only when all pass and all required routes are present."""
    present_routes = {route.route for route in routes}
    missing_routes = [req for req in (required_routes or ()) if req not in present_routes]

    failed_routes = [route.route for route in routes if route.failures]
    overall_status = STATUS_PASS if (not failed_routes and not missing_routes) else STATUS_FAIL

    return {
        "schema": ACCEPTANCE_SCHEMA,
        "status": overall_status,
        "required_routes": list(required_routes or ()),
        "present_routes": list(present_routes),
        "missing_routes": missing_routes,
        "failed_routes": failed_routes,
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
