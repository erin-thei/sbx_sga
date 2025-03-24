import sys

input_reports = snakemake.input
output_report = open(snakemake.output[0], "w")
output_report.write(
    "Sample,Completeness,Contamination,N50,Genome_Size,GC_Content,CDS\n"
)

for file in input_reports:
    file_obj = open(file, "r")
    filelines = file_obj.readlines()
    for line in filelines[1:]:  # Skip header
        line = line.rstrip()
        line_list = line.split("\t")
        sample = line_list[0]
        completeness = line_list[1]
        contamination = line_list[2]
        n50 = line_list[6]
        genome_size = line_list[8]
        gc_content = line_list[9]
        cds = line_list[10]
    output_report.write(
        f"{sample},{completeness},{contamination},{n50},{genome_size},{gc_content},{cds}\n"
    )

output_report.close()
