import pytest
import os
from checkm_f import parse_file, get_stats, write_report


def test_parse_file(tmp_path):
    sample_file = tmp_path / "sample.txt"
    sample_file.write_text(
        """Name	Completeness	Contamination	Completeness_Model_Used	Translation_Table_Used	Coding_Density	Contig_N50	Average_Gene_Length	Genome_Size	GC_Content	Total_Coding_Sequences	Additional_Notes
marc.entero.208	100.0	0.79	Neural Network (Specific Model)	11	0.882	123891	309.07343814646254	5074653	0.51	4834	None
"""
    )
    data_dict = parse_file(sample_file)
    expected_dict = {
        "Name": "marc.entero.208",
        "Completeness": "100.0",
        "Contamination": "0.79",
        "Completeness_Model_Used": "Neural Network (Specific Model)",
        "Translation_Table_Used": "11",
        "Coding_Density": "0.882",
        "Contig_N50": "123891",
        "Average_Gene_Length": "309.07343814646254",
        "Genome_Size": "5074653",
        "GC_Content": "0.51",
        "Total_Coding_Sequences": "4834",
        "Additional_Notes": "None",
    }
    assert data_dict == expected_dict


def test_get_stats():
    expected_dict = {
        "Name": "marc.entero.208",
        "Completeness": "100.0",
        "Contamination": "0.79",
        "Completeness_Model_Used": "Neural Network (Specific Model)",
        "Translation_Table_Used": "11",
        "Coding_Density": "0.882",
        "Contig_N50": "123891",
        "Average_Gene_Length": "309.07343814646254",
        "Genome_Size": "5074653",
        "GC_Content": "0.51",
        "Total_Coding_Sequences": "4834",
        "Additional_Notes": "None",
    }
    assert get_stats(expected_dict) == ("marc.entero.208", "100.0", "0.79")


def test_write_report(tmp_path):
    sample_op = tmp_path / "sample_op.txt"
    sample = "marc.entero.208"
    completeness = "100.0"
    contamination = "0.79"

    op_path = write_report(sample_op, sample, completeness, contamination)
    assert os.path.exists(op_path)

    with open(op_path, "r") as f:
        content = f.read()
        assert f"{sample}\t{completeness}\t{contamination}" in content
