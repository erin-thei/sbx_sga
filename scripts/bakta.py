import sys
from bakta_f import parse_file, get_annotation_stats, write_to_report


report = snakemake.input[0]
output = snakemake.output[0]

filelines = open(report, "r").readlines()
bakta_dict = parse_file(filelines)
(
    crispr_count,
    hypothetical_count,
    trna_count,
    tmrna_count,
    rrna_count,
    cds_count,
    N50,
    genome_size,
    gc
) = get_annotation_stats(bakta_dict)
write_to_report(
    output,
    report,
    crispr_count,
    hypothetical_count,
    trna_count,
    tmrna_count,
    rrna_count,
    cds_count,
    N50,
    genome_size,
    gc
)
