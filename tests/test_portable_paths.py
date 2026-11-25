import shutil
from pathlib import Path

from click.testing import CliRunner

from ultimate.cli import main
from ultimate.config import load_config
from ultimate.paths import project_root


def test_portable_checkout_and_resource_overrides(tmp_path, monkeypatch):
    monkeypatch.delenv("ULTIMATE_ROOT", raising=False)
    assert (project_root() / "config/default.yaml").is_file()
    monkeypatch.chdir(tmp_path)
    checkout = tmp_path / "new checkout with spaces"
    directory = checkout / "config"
    directory.mkdir(parents=True)
    config = directory / "project.yaml"
    config.write_text(
        'project:\n  server_root: ..\n  output_dir: ../results\n'
        'resources:\n  human:\n    gtf: "@resources/GRCh38.gtf"\n'
        '  public_cache: "@project/cache/publicdb"\n'
    )
    resources = tmp_path / "shared references"
    runtime = tmp_path / "site runtime"
    monkeypatch.setenv("ULTIMATE_RESOURCES", str(resources))
    monkeypatch.setenv("ULTIMATE_ROOT", str(runtime))
    loaded = load_config(config).raw
    assert Path(loaded["project"]["server_root"]) == checkout
    assert Path(loaded["project"]["output_dir"]) == checkout / "results"
    assert Path(loaded["resources"]["human"]["gtf"]) == resources / "GRCh38.gtf"
    assert Path(loaded["resources"]["public_cache"]) == checkout / "cache/publicdb"
    config.write_text("{}\n")
    assert Path(load_config(config).raw["project"]["server_root"]) == runtime
    for args, expected in [([], runtime), (["--root", str(checkout)], checkout)]:
        result = CliRunner().invoke(
            main, ["prepare-job", "--config", str(config), "--job-id", "portable_smoke", "--run-mode", "interactive", *args]
        )
        assert result.exit_code == 0, result.output
        assert (expected / "jobs/portable_smoke/config/project.yaml").is_file()


def test_move_config_directory(tmp_path):
    a, b = tmp_path / "A", tmp_path / "B"
    (a / "config").mkdir(parents=True)
    (a / "config/project.yaml").write_text('project:\n  server_root: ..\n  output_dir: ../runs\n')
    shutil.move(a, b)
    result = load_config(b / "config/project.yaml").raw["project"]
    assert Path(result["server_root"]) == b
    assert Path(result["output_dir"]) == b / "runs"
