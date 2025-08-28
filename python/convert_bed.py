import pandas as pd
import sys

def swap_columns(input_file, output_file, header):
    if header:
        df = pd.read_csv(input_file, sep="\t")
    else:
        df = pd.read_csv(input_file, sep="\t", header=None)

    col1, col2 = df.columns[1], df.columns[2]

    # convert col1, col2 to numpy array
    c1, c2 = df[col1].values, df[col2].values

    # swap 
    mask = c1 > c2

    # swap 
    c1_new = np.where(mask, c2, c1) - 1  # left -1
    c2_new = np.where(mask, c1, c2)      # right

    df[col1], df[col2] = c1_new, c2_new

    df_unique = df.drop_duplicates()
    df_unique.to_csv(output_file, sep="\t", index=False, header=False)
    
    print(f"Processing completed. Output file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        swap_columns(input_file, output_file, header)
