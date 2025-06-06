import pytest
import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from amr_f import parse_file, parse_scope, write_long_report

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


def test_parse_file():
    sample_report = os.path.join(test_data_dir, "sample.amrfinder")
    sample, filelines = parse_file(sample_report)
    assert sample == os.path.basename(test_data_dir)


def test_parse_scope():
    sample_report = os.path.join(test_data_dir, "sample.amrfinder")
    files = [sample_report]
    amr_classes = parse_scope(files)
    assert isinstance(amr_classes, dict)


def test_write_long_report(tmp_path):
    sample_dict = {
        "EFFLUX": {"test_data": ["acrF", "emrD"]},
        "BETA-LACTAM": {"test_data": ["blaEC"]},
    }
    output = tmp_path / "output.report"
    report = write_long_report(sample_dict, output)
    with open(report, "r") as op:
        content = op.readlines()
    assert content[0] == "Sample\tGene\tCategory\n"
    assert "test_data\tacrF\tEFFLUX\n" in content
