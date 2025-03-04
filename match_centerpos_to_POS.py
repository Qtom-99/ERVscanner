import pandas as pd
import argparse

# コマンドライン引数の処理
parser = argparse.ArgumentParser(description="Match positions from output with ranges in a BED file.")
parser.add_argument("output_file", help="Path to the output file from the first script.")
parser.add_argument("bed_file", help="Path to the BED file containing ranges.")
parser.add_argument("result_file", help="Path to the result output file.")
args = parser.parse_args()

# 出力ファイルを読み込み
output_data = pd.read_csv(args.output_file, sep=" ", header=None, names=["chr", "position", "group"])

# BEDファイルを読み込み
bed_data = pd.read_csv(args.bed_file, sep="\t", header=None, names=["chr", "start", "end", "id"])

# 各出力行がBEDファイルの範囲内に含まれているか確認
result = []
for _, out_row in output_data.iterrows():
    # 出力行のポジションがBEDファイルの範囲にあるかどうかを確認
    matched_bed = bed_data[
        (bed_data["chr"] == out_row["chr"]) & 
        (bed_data["start"] <= out_row["position"]) & 
        (bed_data["end"] >= out_row["position"])
    ]
    
    # 一致する範囲が見つかった場合、結果リストに追加
    if not matched_bed.empty:
        for _, bed_row in matched_bed.iterrows():
            result.append([out_row["chr"], out_row["position"], out_row["group"], bed_row["id"]])

# 結果を指定されたファイルにタブ区切りで出力
with open(args.result_file, "w") as f:
    for row in result:
        f.write("\t".join(map(str, row)) + "\n")
