"""Slurm entry scripts must find the XieLab runtime layout.

XieLab keeps environments under ``<site>/ultimate_runtime/envs`` and does not
provide conda/mamba/micromamba on ``$PATH``. Every job script therefore has to
discover ``ultimate-core`` the same way instead of assuming
``<checkout>/.conda/envs/ultimate-core`` exists.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SLURM_DIR = REPO_ROOT / "slurm"


def _sbatch_scripts() -> list[Path]:
    return sorted(SLURM_DIR.glob("*.sbatch"))


def test_slurm_dir_has_scripts():
    assert len(_sbatch_scripts()) >= 40


@pytest.mark.parametrize("script", _sbatch_scripts(), ids=lambda path: path.name)
def test_every_sbatch_parses(script: Path):
    """`bash -n` must pass for each job script (no execution)."""
    completed = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr


def test_ultimate_run_discovers_runtime_envs_first():
    text = (SLURM_DIR / "ultimate_run.sbatch").read_text(encoding="utf-8")
    assert 'RUNTIME_ROOT="${ULTIMATE_RUNTIME_ROOT:' in text
    assert '"$ULTIMATE_ENV_ROOT/ultimate-core"' in text
    # The four documented fallbacks appear in order.
    runtime = text.index('"$ULTIMATE_ENV_ROOT/ultimate-core"')
    venv = text.index('"$PROJECT_ROOT/.venv"')
    conda = text.index('"$PROJECT_ROOT/.conda/envs/ultimate-core"')
    assert runtime < venv < conda


def test_ultimate_run_supports_all_activation_modes():
    text = (SLURM_DIR / "ultimate_run.sbatch").read_text(encoding="utf-8")
    assert "ULTIMATE_ENV" in text
    assert "micromamba" in text
    assert "conda.sh" in text or "CONDA_SH" in text
    assert "bin/activate" in text
    assert 'export PATH="$ENV_PREFIX/bin:$PATH"' in text


def test_ultimate_run_exports_environment_markers():
    text = (SLURM_DIR / "ultimate_run.sbatch").read_text(encoding="utf-8")
    assert 'export ULTIMATE_ENV_ROOT=' in text
    assert 'export ULTIMATE_CACHE_ROOT=' in text
    assert 'export ULTIMATE_RUNTIME_ROOT=' in text


def test_ultimate_run_binds_slurm_threads():
    text = (SLURM_DIR / "ultimate_run.sbatch").read_text(encoding="utf-8")
    assert "SLURM_CPUS_PER_TASK" in text
    assert "ULTIMATE_THREADS" in text


def test_no_script_defaults_to_checkout_conda_only():
    offenders = [
        script.name
        for script in _sbatch_scripts()
        if 'CORE_ENV="${ULTIMATE_CORE_ENV:-$PROJECT_ROOT/.conda' in script.read_text(encoding="utf-8")
    ]
    assert offenders == []


def test_core_env_scripts_prefer_runtime_root():
    checked = 0
    for script in _sbatch_scripts():
        text = script.read_text(encoding="utf-8")
        if "CORE_ENV=" not in text:
            continue
        checked += 1
        assert "ULTIMATE_RUNTIME_ROOT" in text, script.name
        # A script may reach the runtime envs via either $RUNTIME_ROOT or
        # $ENV_ROOT, as long as the runtime path is resolved before any
        # in-checkout fallback.
        runtime_markers = (
            '"$RUNTIME_ROOT/envs/ultimate-core"',
            '"$ENV_ROOT/ultimate-core"',
            "$RUNTIME_ROOT/envs/ultimate-core",
            "$ENV_ROOT/ultimate-core",
        )
        runtime = min((text.index(m) for m in runtime_markers if m in text), default=-1)
        assert runtime >= 0, script.name
        conda_markers = ('"$PROJECT_ROOT/.conda/envs/ultimate-core"', "$PROJECT_ROOT/.conda/envs/ultimate-core")
        conda = min((text.index(m) for m in conda_markers if m in text), default=len(text))
        assert runtime < conda, script.name
    assert checked >= 40


def test_setup_envs_script_has_xielab_slurm_directives():
    text = (SLURM_DIR / "xielab_setup_runtime_envs.sbatch").read_text(encoding="utf-8")
    assert "#SBATCH --account=xielab" in text
    assert "#SBATCH --partition=compute_pro" in text
    assert "#SBATCH --qos=normal" in text
    # Building environments must never be auto-submitted from this script: the
    # only permissible mentions are the usage line in comments.
    code_lines = [
        line for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert not [line for line in code_lines if "sbatch " in line]


def test_setup_envs_script_uses_micromamba_only():
    text = (SLURM_DIR / "xielab_setup_runtime_envs.sbatch").read_text(encoding="utf-8")
    assert "micromamba/bin/micromamba" in text
    assert "MAMBA_ROOT_PREFIX" in text
    for key in ("xielab.core.yml", "xielab.scrna.py.yml", "xielab.bio.r.yml"):
        assert key in text


@pytest.mark.parametrize(
    ("filename", "expected_name"),
    [
        ("xielab.core.yml", "ultimate-core"),
        ("xielab.scrna.py.yml", "ultimate-scrna"),
        ("xielab.bio.r.yml", "ultimate-bio-r"),
    ],
)
def test_runtime_env_definitions(filename: str, expected_name: str):
    import yaml

    payload = yaml.safe_load((REPO_ROOT / "envs" / filename).read_text(encoding="utf-8"))
    assert payload["name"] == expected_name
    assert payload["channels"]
    assert payload["dependencies"]


def test_bio_r_env_covers_required_r_packages():
    import yaml

    payload = yaml.safe_load((REPO_ROOT / "envs" / "xielab.bio.r.yml").read_text(encoding="utf-8"))
    deps = " ".join(str(item) for item in payload["dependencies"])
    for package in ("bioconductor-limma", "bioconductor-deseq2", "bioconductor-edger", "r-wgcna", "r-survival"):
        assert package in deps


def test_scrna_env_covers_harmony_and_squidpy():
    import yaml

    payload = yaml.safe_load((REPO_ROOT / "envs" / "xielab.scrna.py.yml").read_text(encoding="utf-8"))
    deps = [str(item) for item in payload["dependencies"]]
    for package in ("harmonypy", "squidpy", "gseapy", "celltypist"):
        assert package in deps


def test_pyproject_ships_runtime_assets():
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert '"share/ultimate/scripts/R" = ["scripts/R/*.R"]' in text
    assert '"share/ultimate/config" = ["config/default.yaml"]' in text
    assert "share/ultimate/slurm" in text
    assert "ultimate_run.sbatch" in text


def test_pyproject_referenced_files_exist():
    """Every data-files entry must point at real files or the wheel build breaks."""
    import tomllib

    payload = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    data_files = payload["tool"]["setuptools"]["data-files"]
    missing: list[str] = []
    for target_dir, patterns in data_files.items():
        for pattern in patterns:
            if not list(REPO_ROOT.glob(pattern)):
                missing.append(f"{target_dir}: {pattern}")
    assert missing == []


def test_setup_runtime_envs_initializes_manifests_and_installs_core():
    text = (REPO_ROOT / "slurm" / "xielab_setup_runtime_envs.sbatch").read_text(encoding="utf-8")
    assert '"$RUNTIME_ROOT/manifests"' in text
    assert 'pip install --no-deps --no-build-isolation -e "$PROJECT_ROOT"' in text
    envs_core = (REPO_ROOT / "envs" / "xielab.core.yml").read_text(encoding="utf-8")
    assert "setuptools" in envs_core
    assert "wheel" in envs_core


def test_e2e_suite_resolves_runtime_cache_samplesheets_and_scrna_env():
    text = (REPO_ROOT / "slurm" / "xielab_e2e_regression_suite.sbatch").read_text(encoding="utf-8")
    assert 'REFERENCE_ROOT="${ULTIMATE_REFERENCE_ROOT:-$RUNTIME_ROOT/cache}"' in text
    assert 'ULTIMATE_E2E_BULK_SAMPLESHEET' in text
    assert 'ULTIMATE_E2E_PROTEOMICS_SAMPLESHEET' in text
    assert '"$SCRN_ENV/bin/python"' in text
