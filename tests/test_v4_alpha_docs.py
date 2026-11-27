"""Docs-boundary checks that must also hold in a standalone checkout.

The independent ``NanqiShen/ultimate`` repository has no parent control repo, so
``AGENTS.md`` is absent there. Resolution is therefore lazy and a missing
control document becomes a skip instead of a collection-time failure.
"""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _agents_root() -> Path | None:
    """Nearest directory holding the control-repo ``AGENTS.md``, if any."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "AGENTS.md").exists():
            return parent
    return None


def _require_agents_root() -> Path:
    root = _agents_root()
    if root is None:
        pytest.skip("no AGENTS.md above the checkout; standalone repository")
    return root


def _ultimate_readme() -> Path:
    agents_root = _agents_root()
    candidates = [REPO_ROOT / "README.md"]
    if agents_root is not None:
        candidates = [agents_root / "README.md", agents_root / "ultimate" / "README.md", *candidates]
    existing = [path for path in candidates if path.exists()]
    for path in existing:
        if "V4 Alpha Customer Delivery Loop" in path.read_text(encoding="utf-8", errors="ignore"):
            return path
    for path in existing:
        if "Ultimate Bioinfo Workbench" in path.read_text(encoding="utf-8", errors="ignore"):
            return path
    pytest.skip("no Ultimate README.md found for docs boundary checks")


def test_agents_documents_v4_customer_delivery_boundary() -> None:
    text = (_require_agents_root() / "AGENTS.md").read_text(encoding="utf-8")

    if "cluster-local operational policy" in text:
        assert "customer-delivery operations require explicit user approval" in text
        assert "production_backend" in text
        assert "customer delivery" in text
        return

    assert "V4 客户交付闭环" in text
    assert "customer_delivery" in text
    assert "internal_rehearsal" in text
    assert "禁止把 `internal_rehearsal` 当作客户交付" in text
    assert "客户版 delivery QA" in text
    assert "batch order scaffold" in text


def test_readme_documents_v4_alpha_customer_delivery_path() -> None:
    text = _ultimate_readme().read_text(encoding="utf-8")

    assert "V4 Alpha Customer Delivery Loop" in text
    assert "ultimate prepare-batch" in text
    assert "v4_alpha_customer_delivery_rehearsal.sbatch" in text
    assert "delivery_mode=customer_delivery_rehearsal" in text
    assert "customer_delivery_rehearsal" in text
    assert "customer_delivery" in text


def test_readme_documents_v4_beta_customer_trial_path() -> None:
    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text(encoding="utf-8")

    assert "V4 Beta Customer Trial" in text
    assert "v4_beta_customer_trial.sbatch" in text
    assert "readme_for_customer.md" in text
    assert "raw upstream evidence is still lightweight" in text
    assert "V4 Beta remains a controlled rehearsal" in text


def test_readme_documents_v4_1_tool_completeness_path() -> None:
    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text(encoding="utf-8")

    assert "V4.1 Tool Completeness" in text
    assert "ultimate tool-completeness" in text
    assert "ultimate order-readiness" in text
    assert "tool_completeness_matrix.tsv" in text
    assert "module_order_readiness_matrix.tsv" in text
    assert "missing_review=false" in text
    assert "sanitization.tsv" in text
    assert "customer_delivery_sanitization.tsv" in text
    assert "customer_package_manifest.tsv" in text


def test_readme_documents_v4_2_raw_to_customer_path() -> None:
    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text(encoding="utf-8")

    assert "V4.2 Raw-to-Customer Trial" in text
    assert "ultimate customer-package" in text
    assert "ultimate batch-status" in text
    assert "v4_2_raw_to_customer_trial.sbatch" in text
    assert "rnaseq_fastq_tiny_counts" in text
    assert "must write a blocked" in text


def test_readme_documents_realistic_pilot_suite() -> None:
    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text(encoding="utf-8")

    assert "Realistic Pilot Suite" in text
    assert "v4_pilot_realistic_suite.sbatch" in text
    assert "v4_pilot_realistic_report.md" in text
    assert "raw-upstream-validated" in text
    assert "delivery-check-ready" in text


def test_docs_boundary_helpers_degrade_without_control_repo(tmp_path, monkeypatch) -> None:
    """A checkout with no parent AGENTS.md must skip, never raise at import time."""
    import importlib.util
    import sys

    import test_v4_alpha_docs as module

    # Standalone checkout: no AGENTS.md anywhere above the copy.
    standalone = tmp_path / "standalone" / "tests"
    standalone.mkdir(parents=True)
    copy = standalone / "test_copy.py"
    copy.write_text(
        Path(module.__file__).read_text(encoding="utf-8").split("def test_")[0],
        encoding="utf-8",
    )
    spec = importlib.util.spec_from_file_location("standalone_alpha_docs", copy)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)  # must not raise

    assert loaded._agents_root() is None
    with pytest.raises(BaseException) as excinfo:
        loaded._require_agents_root()
    assert excinfo.value.__class__.__name__ == "Skipped"

    readme = standalone / "README.md"
    readme.write_text("# Ultimate Bioinfo Workbench\n\nDelivery boundary: rehearsal only.\n", encoding="utf-8")
    # README lookup falls back to the in-repo README when the control repo is missing.
    assert loaded._ultimate_readme().name == "README.md"


def test_standalone_readme_documents_delivery_boundary() -> None:
    """In a standalone checkout, README.md is the delivery-boundary source of truth."""
    if _agents_root() is not None:
        pytest.skip("control repository present; AGENTS.md covers the boundary")
    text = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    lowered = text.lower()
    assert "delivery" in lowered
    assert "rehearsal" in lowered or "customer" in lowered
