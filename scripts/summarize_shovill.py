import os 

input_files = snakemake.input
output_file = snakemake.output[0]
output_obj = open(output_file, "w")
output_obj.write(
    "Sample,Minimum_Contig_Coverage,Avg_Contig_Coverage,Maximum_Contig_Coverage\n"
)

for file in input_files:
    sample = file.split("/")[-1].split(".fa")[0]
    if os.path.getsize(file) > 0: 
        file_obj = open(file, "r")
        filelines = file_obj.readlines()

        total_length = 0
        contig_list = []
        for line in filelines:
            line = line.rstrip()
            if line.startswith(">"):
                line_list = line.split()
                cov = float(line_list[2].split("=")[-1])
                length = float(line_list[1].split("=")[-1])
                contig_cov = cov * length
                contig_list.append(contig_cov)
                total_length += length
            else:
                continue
        total_contig_cov = sum(contig_list)
        min_cov = round(min(contig_list), 2)
        max_cov = round(max(contig_list), 2)
        coverage = round(total_contig_cov / total_length, 2)
        output_obj.write(f"{sample},{min_cov},{coverage},{max_cov}\n")
    else:
        output_obj.write(f"{sample},-,-,-\n")
