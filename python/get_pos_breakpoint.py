import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="各ポジションIDごとに2列目（start）の中央値に15を加えた値をブレークポイントとして出力します。"
    )
    parser.add_argument(
        "input_file",
        help="入力ファイルのパス（例: breakpoints_input.tsv）"
    )
    parser.add_argument(
        "output_file",
        help="出力ファイルのパス（例: breakpoints.tsv）"
    )
    args = parser.parse_args()

    # 入力ファイルはヘッダー無しと仮定し、以下の列名を設定
    # 1列目: 染色体番号 (chr)
    # 2列目: 始点 (start)
    # 3列目: 終点 (end) -- ※ 各行で start との差は常に30
    # 4列目: ポジションID (例: POS1, POS3, …)
    # 5列目: 補助情報（今回は使用しません）
    col_names = ["chr", "start", "end", "posID", "support"]
    df = pd.read_csv(args.input_file, sep="\t", header=None, names=col_names)

    results = []
    # ポジションIDごとにグループ化
    for posID, group in df.groupby("posID"):
        # 染色体番号は同一であると仮定（複数あれば最初の値を使用）
        chrom = group["chr"].iloc[0]
        # start 列の中央値を計算
        median_start = group["start"].median()
        # ブレークポイントは start の中央値に15を加えた値（中央値が小数の場合は四捨五入）
        central_base = int(round(median_start + 15))
        results.append([chrom, central_base, posID])

    # 結果を DataFrame にまとめ、ヘッダー無しのタブ区切りテキストとして出力
    out_df = pd.DataFrame(results, columns=["chr", "central_base", "posID"])
    out_df.to_csv(args.output_file, sep="\t", index=False, header=False)

    print("Breakpoint detection finished. Output file:", args.output_file)

if __name__ == "__main__":
    main()
