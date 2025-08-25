import pandas as pd
import sys

def swap_columns(input_file, output_file, header):
    """
    Read a tab-separated file and swap the values of the second and third columns 
    if the value in the second column is greater than the third column.
    
    :param input_file: Path to the input file
    :param output_file: Path to the output file
    """
    # Load the file
    df = pd.read_csv(input_file, sep='\t', header=header)

    # Swap values if the second column is greater than the third column
    for index, row in df.iterrows():
        col1 = df.columns[1]
        col2 = df.columns[2]

        if row.iloc[1] > row.iloc[2]:
            df.at[index, col1], df.at[index, col2] = row.iloc[2], row.iloc[1]
            row1 = row.iloc[2]  # now in col1
        else:
            row1 = row.iloc[1]  # already the smaller one

        # Subtract 1 from the smaller value (now in col1)
        df.at[index, col1] = row1 - 1

    # Save the modified DataFrame to a new file
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Processing completed. Output file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        swap_columns(input_file, output_file, header)
