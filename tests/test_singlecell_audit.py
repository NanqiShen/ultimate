from __future__ import annotations

from pathlib import Path

import pytest

from ultimate.singlecell_audit import _check_data, run_singlecell_audit


@pytest.fixture(autouse=True)
def isolate_singlecell_audit_env(monkeypatch) -> None:
    monkeypatch.delenv("ULTIMATE_RUNTIME_ROOT", raising=False)
    monkeypatch.delenv("ULTIMATE_ENV_ROOT", raising=False)
    monkeypatch.delenv("ULTIMATE_ROOT", raising=False)


def test_external_data_checks_do_not_use_legacy_paths_by_default(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.delenv("ULTIMATE_SINGLECELL_NSCLC_H5AD", raising=False)
    monkeypatch.delenv("ULTIMATE_SINGLECELL_DNA_BAM_0518", raising=False)
    monkeypatch.delenv("ULTIMATE_SINGLECELL_MTDNA_VARIANTS_0518", raising=False)

    assert _check_data(tmp_path, ("nsclc_h5ad", "dna_bam_0518", "mtdna_0518")) == {
        "nsclc_h5ad": False,
        "dna_bam_0518": False,
        "mtdna_0518": False,
    }


def test_external_data_checks_use_explicit_environment_paths(monkeypatch, tmp_path: Path) -> None:
    nsclc_h5ad = tmp_path / "Lambrechts_all_harmonized.h5ad"
    nsclc_h5ad.write_text("fixture", encoding="utf-8")
    monkeypatch.setenv("ULTIMATE_SINGLECELL_NSCLC_H5AD", str(nsclc_h5ad))

    assert _check_data(tmp_path, ("nsclc_h5ad",)) == {"nsclc_h5ad": True}


def test_singlecell_audit_writes_matrix(tmp_path: Path) -> None:
    root = tmp_path / "ultimate"
    env_root = tmp_path / "ultimate_runtime" / "envs"
    (env_root / "ultimate-core" / "bin").mkdir(parents=True)
    manifest = run_singlecell_audit(root=root, output_dir=tmp_path / "audit")
    assert Path(manifest["capability_matrix"]).exists()
    assert Path(manifest["dependency_report"]).exists()
    assert len(manifest["capabilities"]) >= 13
    assert "delivery_summary" in manifest


def test_genotype_demux_accepts_alternative_command_route(tmp_path: Path) -> None:
    root = tmp_path / "ultimate"
    env_root = tmp_path / "ultimate_runtime" / "envs"
    bin_dir = env_root / "ultimate-genome-mtdna" / "bin"
    bin_dir.mkdir(parents=True)
    for command in ("python", "cellsnp-lite", "vireo"):
        path = bin_dir / command
        path.write_text("#!/usr/bin/env sh\nexit 0\n", encoding="utf-8")
        path.chmod(0o755)
    (env_root / "ultimate-core" / "bin").mkdir(parents=True)
    manifest = run_singlecell_audit(root=root, output_dir=tmp_path / "audit")
    genotype = next(row for row in manifest["capabilities"] if row["capability"] == "genotype_demux")
    assert "souporcell" not in genotype["command_missing"]
