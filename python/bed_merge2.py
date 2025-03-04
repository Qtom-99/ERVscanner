import pandas as pd
import argparse

def merge_bed_regions(input_file_path, output_file_path):
    # BEDファイルの読み込み
    bed_df = pd.read_csv(input_file_path, sep='\t', header=None)

    # 二行ごとにマージ
    merged_rows = []
    for i in range(0, len(bed_df), 2):
        if i+1 < len(bed_df):
            start = min(bed_df.iloc[i, 1], bed_df.iloc[i+1, 1])
            end = max(bed_df.iloc[i, 2], bed_df.iloc[i+1, 2])
            merged_rows.append([bed_df.iloc[i, 0], start, end])
        else:
            merged_rows.append([bed_df.iloc[i, 0], bed_df.iloc[i, 1], bed_df.iloc[i, 2]])

    # 新しいBEDファイルの作成
    merged_bed_df = pd.DataFrame(merged_rows)
    merged_bed_df.to_csv(output_file_path, sep='\t', header=False, index=False)

def main():
    parser = argparse.ArgumentParser(description='Merge BED file regions two rows at a time.')
    parser.add_argument('input_file', help='Path to the input BED file')
    parser.add_argument('output_file', help='Path to the output BED file')

    args = parser.parse_args()

    merge_bed_regions(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
