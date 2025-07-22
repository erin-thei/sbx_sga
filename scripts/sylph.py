import sys
from scripts.sylph_f import parse_file, get_stats, write_report

report = snakemake.input[0]
output = snakemake.output[0]

sylph_stats = parse_file(report)
sample_name, taxo_abundance, contig = get_stats(sylph_stats)
write_report(output, sample_name, taxo_abundance, contig)