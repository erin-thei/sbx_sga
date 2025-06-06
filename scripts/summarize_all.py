from summarize_all_f import summarize_all
import sys

summarize_all(snakemake.input, snakemake.output[0])
