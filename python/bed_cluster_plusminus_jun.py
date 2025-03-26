import sys

def filter_and_save(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Categorize rows based on the value in the 5th column
        rows_by_fifth_column = {}
        for line in infile:
            columns = line.strip().split('\t')
            fifth_column_value = columns[4]

            # Add rows to a dictionary using the 5th column value as the key
            if fifth_column_value in rows_by_fifth_column:
                rows_by_fifth_column[fifth_column_value].append(line)
            else:
                rows_by_fifth_column[fifth_column_value] = [line]

        # Write rows to a new file where rows with the same 5th column value have 4th column values in the order "+" followed by "-"
        for key, rows in rows_by_fifth_column.items():
            plus_minus_values = [row.strip().split('\t')[3] for row in rows]

            # Verify that the order is "+" followed by "-"
            if plus_minus_values == ["+", "-"]:
                for row in rows:
                    outfile.write(row)

if __name__ == "__main__":
    # Parse command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.tsv output_file.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Call a function to extract rows and save them to a new file
    filter_and_save(input_file, output_file)
