import pandas as pd
import sys

def main(file1_path, file2_path, output_path):
    """
    Read two text files and filter rows from the first file 
    based on exact matches in the third column with values from the second file.

    :param file1_path: Path to the first tab-separated file
    :param file2_path: Path to the second newline-separated file
    :param output_path: Path to save the filtered output file
    """
    # Read the first text file (tab-separated)
    df1 = pd.read_csv(file1_path, sep='\t', header=None)

    # Read the second text file (newline-separated)
    with open(file2_path, 'r') as file:
        lines = file.read().splitlines()

    # Convert the second file's content into a set for faster lookup
    search_terms = set(lines)

    # Perform an exact match search on the third column of the first file
    filtered_df = df1[df1[2].isin(search_terms)]

    # Save the filtered results to the output file
    filtered_df.to_csv(output_path, sep='\t', index=False, header=False)
    
    # Print the filtered results
    print(filtered_df)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1> <file2> <output>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_path = sys.argv[3]

    main(file1_path, file2_path, output_path)
