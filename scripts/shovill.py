import sys
from scripts.shovill_f import get_individual_cov, calc_cov_stats, write_to_report

genome = snakemake.input[0]
output = snakemake.output[0]

file_obj = open(genome, "r")
filelines = file_obj.readlines()

if len(filelines) == 0:
    genome_cov = "NA"
    num_contigs = "NA"
    write_to_report(output, genome, genome_cov, num_contigs)
else:
    ctg_stats = get_individual_cov(filelines)
    genome_cov, num_contgs = calc_cov_stats(ctg_stats)
    write_to_report(output, genome, genome_cov, num_contgs)
