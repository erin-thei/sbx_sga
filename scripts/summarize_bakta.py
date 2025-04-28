import os 

input_files = snakemake.input
output_file = snakemake.output[0]
output_obj = open(output_file, "w")
output_obj.write(f"Sample,rRNA,tRNA,tmRNA\n")

for file in input_files:
    sample = file.split("/")[-1].split(".txt")[0]
    if os.path.getsize(file) > 0:
        file_object = open(file, "r")
        filelines = file_object.readlines()
        rrna_count = 0
        trna_count = 0
        tmrna_count = 0
        for line in filelines:
            line = line.rstrip()
            if "rRNA" in line:
                rrna = int(line.split()[-1])
                rrna_count += rrna
            elif "tRNA" in line:
                trna = int(line.split()[-1])
                trna_count += trna
            elif "tmRNA" in line:
                tmrna = int(line.split()[-1])
                tmrna_count += tmrna
            else:
                continue
        output_obj.write(f"{sample},{rrna_count},{trna_count},{tmrna_count}\n")
    else:
        output_obj.write(f"{sample},-,-,-\n")
output_obj.close()
