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
            "modify",
            "-i",
            "-s",
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
            "modify",
            "-i",
            "-s",
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
            "modify",
            "-i",
            "-s",
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

    # Run the test job
    try:
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--profile",
                project_dir,
                "all_sga",
                "--directory",
                temp_dir,
                "-n",
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(log_fp, "logs/")
        shutil.copytree(stats_fp, "stats/")
        sys.exit(e)

    shutil.copytree(log_fp, "logs/")
    shutil.copytree(stats_fp, "stats/")

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam
