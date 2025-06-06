import sys
from scripts.checkm_f import parse_file, get_stats, write_report


report = snakemake.input[0]
output = snakemake.output[0]

checkm_stats = parse_file(report)
sample, completeness, contamination = get_stats(checkm_stats)
write_report(output, sample, completeness, contamination)
