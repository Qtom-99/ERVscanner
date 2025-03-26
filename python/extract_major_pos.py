import argparse
import pandas as pd

def filter_positions(input_file, positions_file, output_file):
    # Read position names separated by newlines
    with open(positions_file, 'r', encoding='utf-8') as f:
        target_positions = set(line.strip() for line in f if line.strip())

    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t', encoding='utf-8')

    # Keep only the target positions for filtering (always retain the sample column)
    filtered_columns = ['Sample'] + [pos for pos in df.columns if pos in target_positions]

    # filtering
    filtered_df = df[filtered_columns]

    # output
    filtered_df.to_csv(output_file, sep='\t', index=False, encoding='utf-8')

def main():
    parser = argparse.ArgumentParser(description="Filter TSV columns based on a list of positions.")
    parser.add_argument('input_file', help="Path to the input TSV file.")
    parser.add_argument('positions_file', help="Path to the file with position names (one per line).")
    parser.add_argument('output_file', help="Path to the output TSV file.")

    args = parser.parse_args()
    filter_positions(args.input_file, args.positions_file, args.output_file)

if __name__ == '__main__':
    main()
