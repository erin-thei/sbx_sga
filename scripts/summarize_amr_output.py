import csv

input_reports = snakemake.input
output_report = snakemake.output[0]

AMR_dict = {}

subclassList = {}
for file in input_reports:
    sample = file.split('/')[-2]
    file_obj = open(file, 'r')
    filelines = file_obj.readlines()
    for line in filelines[1:]:
        line = line.rstrip()
        line_list = line.split('\t')
        scope = line_list[8]
        gene = line_list[5]
                        
        if scope == 'AMR':
            subclass = line_list[11]
            if subclass not in subclassList:
                subclassList[subclass] ={}
            if sample not in subclassList[subclass]:
                subclassList[subclass][sample] = [gene]
            else:
                subclassList[subclass][sample].append(gene)
            if sample in AMR_dict:
                AMR_dict[sample].append(gene)
            else:
                AMR_dict[sample] = [gene]
        else:
            continue

# Extracting all sample names and subclass keys
samples = set()
subclasses = set(subclassList.keys())

for subclass in subclassList:
    samples.update(subclassList[subclass].keys())

# Preparing data structure for TSV
tsv_data = {}
for sample in samples:
    tsv_data[sample] = {subclass: '' for subclass in subclasses}

    for subclass in subclasses:
        if sample in subclassList[subclass]:
            tsv_data[sample][subclass] = '\t'.join(subclassList[subclass][sample])

# Writing to a file
with open(output_report, 'w', newline='') as tsv_file:
    writer = csv.writer(tsv_file, delimiter=',')

    # Writing the header
    header = ['Sample'] + list(subclasses)
    writer.writerow(header)

    # Writing the rows
    for sample, subclass_dict in tsv_data.items():
        row = [sample] + [subclass_dict[subclass] for subclass in subclasses]
        writer.writerow(row)


