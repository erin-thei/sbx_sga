import re 

reports = snakemake.input
output = snakemake.output[0]

mash_annotation_dict = {}

for file in reports:
    sample = file.split('/')[-1].split('_')[0]
    target_species = []
    with open(file,'r') as opened_report:
        mash_lines = opened_report.readlines()
        for report_line in mash_lines[:10]:
            line_list = report_line.rstrip().split('\t')
            median_multiplicity = float(line_list[2])
            if (float(line_list[0]) >= 0.85) and (int(line_list[1].split('/')[0]) >= 100):
                species_line = line_list[-1]
                split_char = re.findall('N[A-Z]_[0-9A-Z]+\.[0-9]',species_line)[0]
                species_split = report_line.split(split_char)[1].lstrip()
                species = ' '.join(species_split.split()[:2])
                if 'phage' in species: # Exclude phages
                    continue
                elif 'Phage' in species: # Exclude phages
                    continue
                elif 'sp.' in species: # Exclude genus-level
                    continue
                else:
                    if median_multiplicity >= 10:
                        target_species.append(species)
    target_set = set(target_species)

    if len(target_set) <= 1:
        mash_annotation_dict[sample] = ''
    else:
        mash_annotation_dict[sample] = 'Contaminated'

with open(output, 'w') as out:
    out.write('Sample,Contamination\n')
    for sample,status in mash_annotation_dict.items():
        if not status == 'Contaminated':
            out.write(f'{sample},None\n')
        else:
            out.write(f'{sample},Contaminated\n')

