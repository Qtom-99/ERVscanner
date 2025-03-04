import argparse
import pandas as pd

def main():
    # 引数の設定
    parser = argparse.ArgumentParser(
        description="サンプルIDを含む TSV ファイルから、全サンプルで0の列および指定された POSITION ID の列を除去します。"
    )
    parser.add_argument(
        "input_file",
        help="入力TSVファイルのパス (例: insertion_data.tsv)"
    )
    parser.add_argument(
        "remove_list_file",
        help="除去すべきPOSITION IDのリストが記載されたテキストファイルのパス (改行区切り)"
    )
    parser.add_argument(
        "output_file",
        help="出力TSVファイルのパス (例: filtered_insertion_data.tsv)"
    )
    args = parser.parse_args()

    # --- 1. TSV ファイルの読み込み ---
    # 1列目が "Sample" で、残りの列が各 POSITION ID の情報と想定
    df = pd.read_csv(args.input_file, sep="\t")
    
    # --- 2. 全サンプルで 0 の列を除去 ---
    # "Sample" 列以外の部分だけを対象にする
    data = df.drop(columns=["Sample"])
    # 各列において、1つでも 0 以外の値があるか確認
    mask = (data != 0).any(axis=0)
    filtered_data = data.loc[:, mask]

    # --- 3. 除去リストに含まれる列の除去 ---
    # remove_list_file 内の各行に "POSx" の形式で列名が記載されているとする
    with open(args.remove_list_file, "r", encoding="utf-8") as f:
        remove_list = [line.strip() for line in f if line.strip()]
    filtered_data = filtered_data.drop(columns=remove_list, errors='ignore')

    # --- 4. 結果の出力 ---
    # サンプルID列とフィルタ済みデータを連結して出力
    result = pd.concat([df["Sample"], filtered_data], axis=1)
    result.to_csv(args.output_file, sep="\t", index=False)

    print("フィルタリングが完了しました。出力ファイル:", args.output_file)

if __name__ == "__main__":
    main()
