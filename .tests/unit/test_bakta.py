import pytest
from bakta_f import parse_file, get_annotation_stats, write_to_report
import os


def test_parse_file(tmp_path):
    sample_file = tmp_path / "sample.txt"
    sample_file.write_text(
        """Sequence(s):
Length: 5074653
Count: 120
GC: 50.8
N50: 123891
N ratio: 0.0
coding density: 88.2
Annotation:
tRNAs: 75
tmRNAs: 1
rRNAs: 5
ncRNAs: 177
ncRNA regions: 52
CRISPR arrays: 0
CDSs: 4743
pseudogenes: 2
hypotheticals: 17
signal peptides: 0
sORFs: 91
gaps: 0
oriCs: 4
oriVs: 0
oriTs: 1
Bakta:
Software: v1.7.0
Database: v5.0
DOI: 10.1099/mgen.0.000685
URL: github.com/oschwengers/bakta
"""
    )
    sample_filelines = open(sample_file, "r").readlines()
    parsed_dict = parse_file(sample_filelines)
    expected_output = {
        "Sequence(s)": "",
        "Length": "5074653",
        "Count": "120",
        "GC": "50.8",
        "N50": "123891",
        "N ratio": "0.0",
        "coding density": "88.2",
        "Annotation": "",
        "tRNAs": "75",
        "tmRNAs": "1",
        "rRNAs": "5",
        "ncRNAs": "177",
        "ncRNA regions": "52",
        "CRISPR arrays": "0",
        "CDSs": "4743",
        "pseudogenes": "2",
        "hypotheticals": "17",
        "signal peptides": "0",
        "sORFs": "91",
        "gaps": "0",
        "oriCs": "4",
        "oriVs": "0",
        "oriTs": "1",
        "Bakta": "",
        "Software": "v1.7.0",
        "Database": "v5.0",
        "DOI": "10.1099/mgen.0.000685",
        "URL": "github.com/oschwengers/bakta",
    }
    assert parsed_dict == expected_output


def test_get_annotation_stats():
    expected_dict = {
        "Sequence(s)": "",
        "Length": "5074653",
        "Count": "120",
        "GC": "50.8",
        "N50": "123891",
        "N ratio": "0.0",
        "coding density": "88.2",
        "tRNAs": "75",
        "tmRNAs": "1",
        "rRNAs": "5",
        "ncRNAs": "177",
        "ncRNA regions": "52",
        "CRISPR arrays": "0",
        "CDSs": "4743",
        "pseudogenes": "2",
        "hypotheticals": "17",
        "signal peptides": "0",
        "sORFs": "91",
        "gaps": "0",
        "oriCs": "4",
        "oriVs": "0",
        "oriTs": "1",
    }
    assert get_annotation_stats(expected_dict) == (
        "0",
        "17",
        "75",
        "1",
        "5",
        "4743",
        "123891",
        "5074653",
        "50.8",
    )


def test_write_to_report(tmp_path):
    output = tmp_path / "test_genome.fasta"
    genome = f"{os.path.splitext(os.path.basename(output))[0]}"
    crispr_count = 0
    hypothetical_count = 17
    trna_count = 75
    tmrna_count = 1
    rrna_count = 5
    cds_count = 4743
    n50 = 123891
    genome_size = 5074653
    gc = 50.8

    op_path = write_to_report(
        output,
        genome,
        crispr_count,
        hypothetical_count,
        trna_count,
        tmrna_count,
        rrna_count,
        cds_count,
        n50,
        genome_size,
        gc,
    )

    assert os.path.exists(op_path)
    with open(op_path, "r") as f:
        content = f.read()
        assert (
            f"{os.path.splitext(os.path.basename(genome))[0]}\t{genome_size}\t{cds_count}\t{n50}\t{rrna_count}\t{trna_count}\t{tmrna_count}\t{crispr_count}\t{hypothetical_count}\t{gc}\n"
            in content
        )
