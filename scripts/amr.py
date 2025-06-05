from scripts.amr_f import parse_scope, write_long_report
import sys

input_reports = snakemake.input
output = snakemake.output[0]

subclass_list = parse_scope(input_reports)
write_long_report(subclass_list, output)
