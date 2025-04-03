import pandas as pd
from functools import reduce

input_files = snakemake.input
output = snakemake.output[0]


def process_filelines(file, tool, master_list):

    if tool == "mlst":
        df = pd.read_csv(file, header=None, sep="\t")
        df.columns = [
            "Sample",
            "Schema",
            "ST",
            "Allele_1",
            "Allele_2",
            "Allele_3",
            "Allele_4",
            "Allele_5",
            "Allele_6",
            "Allele_7",
        ]
        df["Sample"] = df["Sample"].apply(lambda x: x.rstrip(".fa"))
        master_list.append(df)

    else:
        df = pd.read_csv(file)
        master_list.append(df)


master_list = []
for file in input_files:
    tool = file.split("/")[-1].split(".report")[0]
    process_filelines(file, tool, master_list)

final_df = reduce(
    lambda left, right: pd.merge(left, right, on="Sample", how="outer"), master_list
)
final_df.to_csv(output, index=False)
