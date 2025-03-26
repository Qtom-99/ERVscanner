import pandas as pd
import sys

def swap_columns(input_file, output_file):
    """
    Read a tab-separated file and swap the values of the second and third columns 
    if the value in the second column is greater than the third column.
    
    :param input_file: Path to the input file
    :param output_file: Path to the output file
    """
    # Load the file
    df = pd.read_csv(input_file, sep='\t')

    # Swap values if the second column is greater than the third column
    for index, row in df.iterrows():
        if row[1] > row[2]:
            df.at[index, df.columns[1]], df.at[index, df.columns[2]] = row[2], row[1]

    # Save the modified DataFrame to a new file
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Processing completed. Output file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        swap_columns(input_file, output_file)
