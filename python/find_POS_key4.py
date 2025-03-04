import sys
import os
import pandas as pd
from multiprocessing import Pool, cpu_count

# 引数の数が正しいことを確認
if len(sys.argv) != 5:
    print("Usage: python script.py file1.txt file2.txt output.txt num_cores")
    sys.exit(1)

file1_path = sys.argv[1]
file2_path = sys.argv[2]
output_path = sys.argv[3]
num_cores = int(sys.argv[4])

# 使用可能なCPUコア数を取得
max_cores = cpu_count()

# 指定されたコア数が使用可能な範囲内であることを確認
if num_cores < 1 or num_cores > max_cores:
    print(f"Invalid number of cores specified. Please choose a number between 1 and {max_cores}.")
    sys.exit(1)

# ファイル１の名前を取得
file1_name = os.path.basename(file1_path)

# ファイル１とファイル２の読み込み
df1 = pd.read_csv(file1_path, sep="\t", header=None)
df2 = pd.read_csv(file2_path, sep="\t", header=None)

# ファイル2のデータを辞書に格納
data_dict = {}
for _, row in df2.iterrows():
    chrom = row[1]
    if chrom not in data_dict:
        data_dict[chrom] = []
    data_dict[chrom].append((row[2], row[3], row))

# マッチング関数
def match_rows(row1):
    chrom_file1 = row1[1]
    pos_file1 = row1[2]
    
    if chrom_file1 in data_dict:
        for start_pos, end_pos, data in data_dict[chrom_file1]:
            if start_pos <= pos_file1 <= end_pos:
                return list(data[:4]) + list(row1) + [file1_name]
    
    return ["uniq", "uniq", "uniq", "uniq"] + list(row1) + [file1_name]

# 並列処理
if __name__ == "__main__":
    with Pool(processes=num_cores) as pool:
        result = pool.map(match_rows, [row for _, row in df1.iterrows()])
    
    output_df = pd.DataFrame(result)
    output_df.to_csv(output_path, sep="\t", index=False, header=False)
