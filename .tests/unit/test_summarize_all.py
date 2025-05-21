import pandas
from summarize_all_f import process_filelines, summarize_all


def test_process_filelines(tmp_path):
    # Create a temporary file to simulate the input file
    input_file = tmp_path / "test_input.txt"
    with open(input_file, "w") as f:
        f.write("Sample\tSchema\tST\tAllele_1\tAllele_2\n")
        f.write("sample1\tschema1\tST1\tallele1\tallele2\n")
        f.write("sample2\tschema2\tST2\tallele3\tallele4\n")

    # Create a master list to store the DataFrames
    master_list = []

    # Call the function with the test input
    process_filelines(input_file, "mlst", master_list)

    # Check if the DataFrame is created correctly
    assert len(master_list) == 1
    assert isinstance(master_list[0], pandas.DataFrame)
    assert master_list[0].shape == (3, 10)  # 2 rows, 10 columns
