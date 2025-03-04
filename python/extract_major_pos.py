import argparse
import pandas as pd

def filter_positions(input_file, positions_file, output_file):
    # ポジション名を改行区切りで読み込む
    with open(positions_file, 'r', encoding='utf-8') as f:
        target_positions = set(line.strip() for line in f if line.strip())

    # 入力TSVファイルを読み込む
    df = pd.read_csv(input_file, sep='\t', encoding='utf-8')

    # フィルター対象ポジションを保持（サンプル列は常に残す）
    filtered_columns = ['Sample'] + [pos for pos in df.columns if pos in target_positions]

    # フィルタリング
    filtered_df = df[filtered_columns]

    # 結果を出力
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