import os
import pytest
import shutil
import subprocess as sp
import sys
import tempfile
from pathlib import Path


@pytest.fixture
def setup():
    temp_dir = Path(tempfile.mkdtemp())

    reads_fp = Path(".tests/data/reads/").resolve()

    project_dir = temp_dir / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_sga: {{mash_ref: '{temp_dir}/dummy.msh'}}"
    Path(temp_dir / "dummy.msh").touch()

    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_sga: {{checkm_ref: '{temp_dir}/dummy.1.dmnd'}}"
    Path(temp_dir / "dummy.1.dmnd").touch()

    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"sbx_sga: {{bakta_ref: '{temp_dir}/bakta/db/'}}"
    Path(temp_dir / "bakta/db/").mkdir(parents=True, exist_ok=True)

    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    os.environ["SUNBEAM_EXTENSIONS"] = str(Path("extensions/").resolve())

    # DEBUG
    from sunbeam import EXTENSIONS_DIR

    print("EXTENSIONS_DIR: ", EXTENSIONS_DIR())
    for ext in EXTENSIONS_DIR().iterdir():
        print("EXTENSIONS_DIR: ", ext)
        if ext.is_dir():
            for file in ext.iterdir():
                print("EXTENSIONS_DIR: ", file)

    print(Path("extensions/").resolve())
    for ext in Path("extensions/").iterdir():
        print("WRONG: ", ext)
        if ext.is_dir():
            for file in ext.iterdir():
                print("WRONG: ", file)

    sbx_proc = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_sga",
            "--directory",
            temp_dir,
            "-n",
        ],
        capture_output=True,
        text=True,
    )

    print("STDOUT: ", sbx_proc.stdout)
    print("STDERR: ", sbx_proc.stderr)

    try:
        shutil.copytree(log_fp, "logs/")
        shutil.copytree(stats_fp, "stats/")
    except FileNotFoundError:
        print("No logs or stats directory found.")
        Path("logs/").mkdir(parents=True, exist_ok=True)
        Path("stats/").mkdir(parents=True, exist_ok=True)
        Path("logs/file").touch()
        Path("stats/file").touch()

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp, sbx_proc


def test_dry_run(run_sunbeam):
    output_fp, benchmarks_fp, proc = run_sunbeam

    assert proc.returncode == 0, f"Sunbeam run failed with error: {proc.stderr}"
