import sys
import pandas as pd

def filter_file(input_file_path, output_file_path, min_count):
    # ファイルの読み込み
    df = pd.read_csv(input_file_path, sep='\t', header=None)

    # 7列目の値を基準に行をフィルタリング
    filtered_df = df[df.groupby(6)[6].transform('count') >= min_count]

    # 結果をファイルに書き込み
    filtered_df.to_csv(output_file_path, sep='\t', index=False, header=False)

if __name__ == "__main__":
    # コマンドライン引数の取得
    if len(sys.argv) != 4:
        print("Usage: python bed_curate.py input_file.tsv output_file.tsv min_count")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    try:
        min_count = int(sys.argv[3])
    except ValueError:
        print("min_count must be an integer.")
        sys.exit(1)

    # ファイルのフィルタリング
    filter_file(input_file, output_file, min_count)
