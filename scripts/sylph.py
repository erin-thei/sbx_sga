import sys

sys.stderr.write("Starting Sylph processing...\n")
from scripts.sylph_f import parse_file, get_stats, write_report

report = snakemake.input[0]
output = snakemake.output[0]

sys.stderr.write(f"Processing report: {report}\n")
sylph_stats = parse_file(report)
sys.stderr.write(f"Parsed data: {sylph_stats}\n")
sample_name, taxo_abundance, contig = get_stats(sylph_stats)
sys.stderr.write("Obtained stats!\n")
write_report(output, sample_name, taxo_abundance, contig)
sys.stderr.write(f"Report written to: {output}\n")
