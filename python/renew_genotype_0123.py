import pandas as pd
import os
import sys
import glob

def update_table(tsv_file, sample_files_dir, output_file):
    # 1. TSVファイルの読み込み（行：サンプルID、列：POSID）
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    
    # 2. 指定ディレクトリ内のジェノタイピングファイルを取得
    genotype_files = glob.glob(os.path.join(sample_files_dir, "*_n_of_covreads_noclip.txt"))
    print(f"Found {len(genotype_files)} genotype files in {sample_files_dir}")

    # 3. 各サンプルのファイルを処理
    for file_path in genotype_files:
        base_name = os.path.basename(file_path)
        # ファイル名は "SAMPLEID_n_of_covreads_noclip.txt" であることを想定
        if not base_name.endswith("_n_of_covreads_noclip.txt"):
            continue
        sample_id = base_name.split("_n_of_covreads_noclip.txt")[0]
        print(f"Processing sample: {sample_id}")

        if sample_id not in df.index:
            print(f"Warning: Sample ID {sample_id} not found in TSV. Skipping this sample.")
            continue

        # 対象サンプルの元の行（各POSIDの元の値）を退避
        original_row = df.loc[sample_id].copy()
        # このサンプルファイル内で情報が得られたPOSIDの集合
        observed_posids = set()

        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                # 少なくとも5カラム（4列目：POSID、5列目：score）が必要
                if len(parts) < 5:
                    continue  # 4列目にPOSIDがない行は情報が得られないのでスキップ
                pos_id = parts[3]  # 4列目（POSID）
                try:
                    score = float(parts[4])
                except ValueError:
                    continue

                # この行で対象としたPOSIDを記録
                observed_posids.add(pos_id)
                # TSVの列にそのPOSIDが存在する場合のみ更新
                if pos_id not in df.columns:
                    continue

                # ジェノタイピング結果に基づいて更新
                #  score <= 2  → ホモ：更新値は 2
                #  score >= 3  → ヘテロ：更新値は 1
                if score <= 2:
                    df.loc[sample_id, pos_id] = 2
                elif score >= 3:
                    df.loc[sample_id, pos_id] = 1

        # 4. 【新ルール】：
        # もし、ジェノタイピングファイル内にそのPOSIDの情報がなかった場合で、
        # 元のTSV上でそのPOSIDが 1（挿入あり）となっていたなら、更新値を 3 にする
        for pos in df.columns:
            if pos not in observed_posids and original_row.get(pos, 0) == 1:
                df.loc[sample_id, pos] = 3

    # 5. 更新後のTSVを保存
    df.to_csv(output_file, sep="\t")
    print(f"Updated table saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python update_table.py <tsv_file> <sample_files_dir> <output_file>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    sample_files_dir = sys.argv[2]
    output_file = sys.argv[3]

    update_table(tsv_file, sample_files_dir, output_file)
