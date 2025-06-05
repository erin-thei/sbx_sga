import pandas
from summarize_all_f import process_filelines, summarize_all


def test_process_filelines(tmp_path):
    # Create a temporary file to simulate the input file
    input_file = tmp_path / "test_input.txt"
    with open(input_file, "w") as f:
        f.write("sample1.fa\tsaureus\t398\t1 2 3 4 5 6 7\n")

    # Create a master list to store the DataFrames
    master_list = []

    # Call the function with the test input
    process_filelines(input_file, "mlst", master_list)

    # Check if the DataFrame is created correctly
    assert len(master_list) == 1
    assert isinstance(master_list[0], pandas.DataFrame)
    assert master_list[0].shape == (1, 4)  # 1 rows, 4 columns
