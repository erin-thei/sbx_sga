import pandas as pd
from functools import reduce


def process_filelines(fp, tool, master_list):
    if tool == "mlst":
        colnames = [
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
        df = pd.read_csv(fp, sep="\t", names=colnames, header=None, on_bad_lines="warn")
        df["Sample"] = df["Sample"].apply(lambda x: x.rstrip(".fa"))
        master_list.append(df)

    else:
        df = pd.read_csv(fp)
        master_list.append(df)


def summarize_all(input_files, output):
    if not input_files:
        # Handle empty input_files by writing an empty CSV with just the header
        empty_df = pd.DataFrame(columns=["Sample"])
        empty_df.to_csv(output, index=False)
        return
    master_list = []
    for fp in input_files:
        tool = fp.split("/")[-1].split(".report")[0]
        process_filelines(fp, tool, master_list)

    final_df = reduce(
        lambda left, right: pd.merge(left, right, on="Sample", how="outer"), master_list
    )
    final_df.to_csv(output, index=False)
