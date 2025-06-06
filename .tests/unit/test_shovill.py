import pytest
import os
from shovill_f import get_individual_cov, calc_cov_stats, write_to_report


def test_get_individual_cov():
    example_fasta = [
        ">contig00001 len=129214 cov=78.8 corr=0 origname=Contig_17_78.768_pilon sw=shovill-skesa/1.1.0 date=20250514",
        "AGAAGAGAAAAAAGAAGACGAAGATAAATAATTATCTTTAAAAGTTAATGTATCAAATGA",
        ">contig00039 len=24105 cov=86.0 corr=0 origname=Contig_22_85.9577_pilon sw=shovill-skesa/1.1.0 date=20250514",
        "AGAAGAGAAAAAAGAAGACGAAGATAAATAATTATCTTTAAAAGTTAATGTATCAAATGA",
    ]
    stats = get_individual_cov(example_fasta)
    expected = [[1, 78.8, 129214], [2, 86.0, 24105]]
    assert stats == expected


def test_calc_cov_stats():
    all_ctgs = [[">contig00001", 78.8, 129214], [">contig00039", 86.0, 24105]]
    genome_cov, ctg_count = calc_cov_stats(all_ctgs)
    assert genome_cov == 79.93
    assert ctg_count == 2


def test_write_to_report(tmp_path):
    genome_cov = 79.3
    ctg_count = 2
    sample_path = tmp_path / "fasta.fa"
    sample_path.write_text(
        f">contig00001 len=129214 cov=78.8 corr=0 origname=Contig_17_78.768_pilon sw=shovill-skesa/1.1.0 date=20250514"
    )

    op_path = tmp_path / "shovill.test"
    op_report = write_to_report(op_path, str(sample_path), genome_cov, ctg_count)
    assert os.path.exists(op_report)

    with open(op_report) as f:
        content = f.read().strip()
    expected_sample_name = "fasta"

    assert content == f"{expected_sample_name}\t{genome_cov}\t{ctg_count}"
