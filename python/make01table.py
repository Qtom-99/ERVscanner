import pandas as pd
import argparse

def create_position_sample_matrix(data_file, samples_file, output_file):
    # データを読み込む
    df = pd.read_csv(data_file, sep='\s+', header=None, names=['Position', 'Sample'])

    # サンプルリストを読み込む
    with open(samples_file, 'r') as f:
        samples = f.read().splitlines()

    # Position列をユニーク化し、自然順にソートする (e.g., POS1, POS2, ..., POS10)
    df['Position'] = pd.Categorical(df['Position'], ordered=True, categories=sorted(df['Position'].unique(), key=lambda x: int(x[3:])))

    # Positionを再マッピング (e.g., POS1, POS2, ..., POSn)
    position_mapping = {old_pos: f'POS{i+1}' for i, old_pos in enumerate(df['Position'].cat.categories)}
    df['Position'] = df['Position'].map(position_mapping)



    # サンプルリストをデータフレームに変換して全サンプルを含める
    samples_df = pd.DataFrame(samples, columns=['Sample'])

    # ピボットテーブルを作成 (observed=Falseを明示的に指定して警告を抑制)
    pivot_df = df.pivot_table(index='Sample', columns='Position', aggfunc='size', fill_value=0, observed=False)
    pivot_df = (pivot_df > 0).astype(int)  # データを0/1に変換

    # 全サンプルを含むようにマージ (欠損値は0で補完)
    pivot_df = samples_df.merge(pivot_df, on='Sample', how='left').fillna(0).set_index('Sample')
    pivot_df = pivot_df.astype(int)  # NaNをint型の0に変換

    # タブ区切りのCSVとして保存
    pivot_df.to_csv(output_file, sep='\t')

    return pivot_df

def main():
    parser = argparse.ArgumentParser(description="Create a position-sample matrix")
    parser.add_argument('data_file', help='Path to the data file')
    parser.add_argument('samples_file', help='Path to the samples file')
    parser.add_argument('output_file', help='Path to the output file')
    args = parser.parse_args()

    create_position_sample_matrix(args.data_file, args.samples_file, args.output_file)

if __name__ == "__main__":
    main()
