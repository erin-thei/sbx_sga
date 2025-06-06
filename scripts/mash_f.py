import re
import os


def open_report(report):
    sample_name = os.path.basename(report.split("_sorted_winning.tab")[0])
    with open(report, "r") as report_obj:
        filelines = report_obj.readlines()
    top_lines = filelines[:20]
    return sample_name, top_lines


def process_mash_line(line):
    line_list = line.rstrip().split("\t")
    species_line = line_list[-1]
    split_char = re.findall("N[A-Z]_[0-9A-Z]+\.[0-9]", species_line)[0]
    species_split = line.split(split_char)[1].lstrip()
    species = " ".join(species_split.split()[:2])
    median_multiplicity = float(line_list[2])
    identity = float(line_list[0])
    hits = int(line_list[1].split("/")[0])
    return species, median_multiplicity, identity, hits


def parse_report(top_lines):
    target_species = []

    # Assumption is that the first hit in the sorted report does not contain a phage.
    top_species, top_median_multiplicity, top_identity, top_hits = process_mash_line(
        top_lines[0]
    )
    if (top_identity >= 0.85) and (top_hits >= 100):
        target_species.append(top_species)
    else:
        pass

    # Set the threshold for median multiplicity to 5% of the top hit.
    threshold = 0.05 * top_median_multiplicity

    # Iterate through the rest of the hits, excluding phages and genus-level hits.
    # Also, exclude hits with identity < 0.85 and hits < 100.
    for line in top_lines[1:]:
        species, median_multiplicity, identity, hits = process_mash_line(line)
        if (identity >= 0.85) and (hits >= 100):
            if "phage" in species:  # Exclude phages
                continue
            elif "Phage" in species:  # Exclude phages
                continue
            elif "sp." in species:  # Exclude genus-level
                continue
            else:
                if median_multiplicity >= threshold:
                    target_species.append(species)
    target_set = set(target_species)
    return target_set


def contamination_call(target_set):
    mash_dict = {}
    if len(target_set) <= 1:
        mash_dict["NA"] = ""
    else:
        species = " ".join(sorted(target_set))
        mash_dict["Contaminated"] = species
    return mash_dict


def write_report(output, sample_name, mash_dict):
    # Expecting that the dictionary is just one key-value pair, so need to check that
    if len(mash_dict) == 1:
        status = list(mash_dict.keys())[0]
    else:
        # Raise error if dictionary is not proper length
        pass
    with open(output, "w") as out:
        if status == "Contaminated":
            contaminated_spp = mash_dict[status]
            out.write(f"{sample_name}\tContaminated\t{contaminated_spp}\n")
        else:
            out.write(f"{sample_name}\tNA\tNA\n")
    return output
