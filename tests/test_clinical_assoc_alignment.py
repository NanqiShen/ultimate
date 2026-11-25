"""clinical_assoc alignment and production-mode input guards.

Two regressions are pinned here:

1. ``_clinical_assoc_tables`` used to align clinical values by row *position*
   (``reindex(range(len(values)))``), so any clinical table whose ``sample_id``
   order differed from the matrix columns produced correlations between
   mismatched samples. Alignment must key on ``sample_id``.
2. A ``run_mode=production`` job with a missing input matrix silently fell back
   to a randomly generated demo matrix, producing delivery-looking numbers with
   no evidence. Production must block instead.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ultimate.bulk import _clinical_assoc_tables, require_real_inputs_for_production
from ultimate.config import normalize_config
from ultimate.modules.runner import run_module


def _write_clinical(path: Path, sample_ids, values, extra=None) -> pd.DataFrame:
    frame = pd.DataFrame({"sample_id": list(sample_ids), "score": list(values)})
    if extra:
        for key, column in extra.items():
            frame[key] = column
    frame.to_csv(path, sep="\t", index=False)
    return frame


def _read_assoc(tables_dir: Path) -> pd.DataFrame:
    return pd.read_csv(tables_dir / "clinical_feature_associations.tsv", sep="\t")


def test_clinical_assoc_aligns_shuffled_sample_ids(tmp_path: Path) -> None:
    sample_ids = ["S1", "S2", "S3", "S4", "S5"]
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    matrix = pd.DataFrame(
        [values],
        index=["CLIN_0001"],
        columns=sample_ids,
    )
    clinical = pd.DataFrame({"sample_id": ["S3", "S1", "S5", "S2", "S4"], "score": [3.0, 1.0, 5.0, 2.0, 4.0]})

    _clinical_assoc_tables(matrix, clinical, tmp_path)
    assoc = _read_assoc(tmp_path)
    row = assoc[assoc["clinical_variable"] == "score"].iloc[0]

    # A perfect match can only appear when alignment happens on sample_id.
    assert row["correlation"] == pytest.approx(1.0)
    assert int(row["n_samples"]) == 5


def test_clinical_assoc_detects_misalignment_that_positional_index_hides(tmp_path: Path) -> None:
    """A reversed clinical order must give -1, not +1 as positional alignment would."""
    sample_ids = ["S1", "S2", "S3", "S4"]
    matrix = pd.DataFrame([[1.0, 2.0, 3.0, 4.0]], index=["CLIN_0002"], columns=sample_ids)
    clinical = pd.DataFrame({"sample_id": ["S1", "S2", "S3", "S4"], "score": [4.0, 3.0, 2.0, 1.0]})

    _clinical_assoc_tables(matrix, clinical, tmp_path)
    row = _read_assoc(tmp_path).iloc[0]

    assert row["correlation"] == pytest.approx(-1.0)


def test_clinical_assoc_tolerates_missing_and_extra_samples(tmp_path: Path) -> None:
    sample_ids = ["S1", "S2", "S3", "S4", "S5"]
    matrix = pd.DataFrame([[5.0, 4.0, 3.0, 2.0, 1.0]], index=["CLIN_0003"], columns=sample_ids)
    # S3 is absent here, and S9 is unknown to the matrix.
    clinical = pd.DataFrame(
        {"sample_id": ["S1", "S2", "S4", "S5", "S9"], "score": [5.0, 4.0, 2.0, 1.0, 99.0]}
    )

    _clinical_assoc_tables(matrix, clinical, tmp_path)
    row = _read_assoc(tmp_path).iloc[0]

    assert row["correlation"] == pytest.approx(1.0)
    assert int(row["n_samples"]) == 4


def test_clinical_assoc_skips_constant_and_sparse_variables(tmp_path: Path) -> None:
    sample_ids = ["S1", "S2", "S3", "S4"]
    matrix = pd.DataFrame([[1.0, 2.0, 3.0, 4.0]], index=["CLIN_0004"], columns=sample_ids)
    clinical = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "constant": [7.0, 7.0, 7.0, 7.0],
            "sparse": [1.0, np.nan, np.nan, np.nan],
        }
    )

    _clinical_assoc_tables(matrix, clinical, tmp_path)
    assoc = _read_assoc(tmp_path)

    assert "constant" not in set(assoc["clinical_variable"])
    assert "sparse" not in set(assoc["clinical_variable"])


def test_clinical_assoc_records_no_variable_row_when_nothing_usable(tmp_path: Path) -> None:
    matrix = pd.DataFrame([[1.0, 2.0]], index=["CLIN_0005"], columns=["S1", "S2"])
    clinical = pd.DataFrame({"sample_id": ["S1", "S2"], "score": [3.0, 3.0]})

    _clinical_assoc_tables(matrix, clinical, tmp_path)
    assoc = _read_assoc(tmp_path)

    assert list(assoc["clinical_variable"]) == ["none"]
    assert int(assoc.iloc[0]["n_samples"]) == 0


# --- production input guard --------------------------------------------------


def _production_config(module_name: str, tmp_path: Path, **module_cfg) -> dict:
    config = {
        "project": {"name": "guard", "organism": "human", "run_mode": "production", "output_dir": str(tmp_path / "runs")},
        "design": {"condition_column": "condition", "control": "control", "case": "treated"},
        "modules": {module_name: module_cfg},
        "samples": [],
    }
    return normalize_config(config, tmp_path)


def test_production_missing_matrix_blocks(tmp_path: Path) -> None:
    config = _production_config("rnaseq", tmp_path, input_matrix=str(tmp_path / "missing.tsv"))

    with pytest.raises(ValueError, match="refusing to generate a demo matrix"):
        require_real_inputs_for_production(config, "rnaseq", config["modules"]["rnaseq"], config["modules"]["rnaseq"].get("input_matrix"))


def test_production_without_any_matrix_declaration_blocks(tmp_path: Path) -> None:
    config = _production_config("rnaseq", tmp_path)

    with pytest.raises(ValueError, match="no matrix path declared"):
        require_real_inputs_for_production(config, "rnaseq", config["modules"]["rnaseq"], None)


def test_production_accepts_existing_matrix(tmp_path: Path) -> None:
    matrix = tmp_path / "counts.tsv"
    matrix.write_text("feature_id\tS1\tS2\nG1\t1\t2\n", encoding="utf-8")
    config = _production_config("rnaseq", tmp_path, input_matrix=str(matrix))

    require_real_inputs_for_production(config, "rnaseq", config["modules"]["rnaseq"], str(matrix))


def test_production_accepts_raw_output_matrix_alias(tmp_path: Path) -> None:
    matrix = tmp_path / "upstream.tsv"
    matrix.write_text("feature_id\tS1\tS2\nG1\t1\t2\n", encoding="utf-8")
    config = _production_config("rnaseq", tmp_path, raw={"output_matrix": str(matrix)})

    require_real_inputs_for_production(config, "rnaseq", config["modules"]["rnaseq"], None)


def test_production_explicit_demo_is_allowed(tmp_path: Path) -> None:
    config = _production_config("rnaseq", tmp_path, is_demo=True)

    require_real_inputs_for_production(config, "rnaseq", config["modules"]["rnaseq"], None)


def test_interactive_missing_matrix_is_allowed(tmp_path: Path) -> None:
    config = _production_config("rnaseq", tmp_path, input_matrix=str(tmp_path / "missing.tsv"))
    config["project"]["run_mode"] = "interactive"

    require_real_inputs_for_production(config, "rnaseq", config["modules"]["rnaseq"], config["modules"]["rnaseq"].get("input_matrix"))


def test_runner_blocks_production_run_without_input(tmp_path: Path) -> None:
    """The module runner must not silently fabricate a matrix in production."""
    config = _production_config("rnaseq", tmp_path)
    samples = pd.DataFrame(
        {
            "sample_id": ["C1", "C2", "T1", "T2"],
            "condition": ["control", "control", "treated", "treated"],
        }
    )

    with pytest.raises(ValueError, match="refusing to generate a demo matrix"):
        run_module(module_name="rnaseq", config=config, output_dir=tmp_path / "runs" / "guard", samples=samples)


def test_runner_interactive_run_still_uses_demo_matrix(tmp_path: Path) -> None:
    """Interactive smoke runs keep working, so the guard does not break demos."""
    config = _production_config("clinical_assoc", tmp_path)
    config["project"]["run_mode"] = "interactive"
    samples = pd.DataFrame(
        {
            "sample_id": ["C1", "C2", "T1", "T2"],
            "condition": ["control", "control", "treated", "treated"],
        }
    )

    manifest = run_module(module_name="clinical_assoc", config=config, output_dir=tmp_path / "runs" / "demo", samples=samples)

    assert manifest["module"] == "clinical_assoc"
