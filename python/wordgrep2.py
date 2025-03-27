import pandas as pd
import sys
import csv

def main(file1_path, file2_path, output_path):
    """
    Read two text files and filter rows from the first file 
    based on exact matches in the third column with values from the second file.

    :param file1_path: Path to the first tab-separated file
    :param file2_path: Path to the second newline-separated file
    :param output_path: Path to save the filtered output file
    """
    with open(file2_path, 'r') as file:
        search_terms = set(file.read().splitlines())

    # Open the input and output files
    with gzip.open(file1_path, 'rt') as infile, open(output_path, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Process line by line to minimize memory usage
        for row in reader:
            if len(row) > 2 and row[2] in search_terms:
                writer.writerow(row)
    
    # Print the filtered results
    print("File processed successfully")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1> <file2> <output>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_path = sys.argv[3]

    main(file1_path, file2_path, output_path)
