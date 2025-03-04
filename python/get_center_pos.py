import pandas as pd
import argparse

# コマンドライン引数の処理
parser = argparse.ArgumentParser(description="Extract minimum values from BED file based on matching group.")
parser.add_argument("input_file", help="Path to the input BED file.")
parser.add_argument("output_file", help="Path to the output file.")
args = parser.parse_args()

# BEDファイルを読み込み
data = pd.read_csv(args.input_file, sep="\t", header=None, names=["chr", "start", "end", "strand", "group"])

# グループごとに+と-の行を取得
plus_strands = data[data["strand"] == "+"]
minus_strands = data[data["strand"] == "-"]

# 最小値を計算して結果リストに保存
result = []
for group in data["group"].unique():
    plus_row = plus_strands[plus_strands["group"] == group]
    minus_row = minus_strands[minus_strands["group"] == group]
    
    if not plus_row.empty and not minus_row.empty:
        # +側のend列の値と-側のstart列の値を取得し、最小値を計算
        plus_value = int(plus_row.iloc[0]["end"])
        minus_value = int(minus_row.iloc[0]["start"])
        min_value = min(plus_value, minus_value)
        result.append([plus_row.iloc[0]["chr"], min_value, group])

# 結果を指定されたファイルに出力
with open(args.output_file, "w") as f:
    for row in result:
        f.write(f"{row[0]} {row[1]} {row[2]}\n")
