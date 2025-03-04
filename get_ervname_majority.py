#!/usr/bin/env python3
import argparse
from collections import defaultdict, Counter

def main():
    parser = argparse.ArgumentParser(
        description="TSVファイル（POSID, ERVNAME, ERVCLASS）から、各POSIDごとに最も頻出するERVNAMEを出力します。"
    )
    parser.add_argument("input_file", help="入力TSVファイルのパス")
    parser.add_argument("output_file", help="出力ファイルのパス")
    args = parser.parse_args()

    # POSIDごとにERVNAMEの出現回数を数えるための辞書（Counter を使用）
    posid_counters = defaultdict(Counter)

    # 入力ファイルの読み込み（タブ区切り）
    with open(args.input_file, "r", encoding="utf-8") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            # 少なくともPOSIDとERVNAMEの2列があることを確認
            if len(fields) < 2:
                continue
            posid = fields[0]
            ervname = fields[1]
            posid_counters[posid][ervname] += 1

    # 各POSIDごとに最も頻出するERVNAMEを決定し、出力ファイルに書き出す
    with open(args.output_file, "w", encoding="utf-8") as outfile:
        # 出力フォーマット：POSID<TAB>最も多いERVNAME
        for posid, counter in posid_counters.items():
            if counter:
                most_common_name = counter.most_common(1)[0][0]
            else:
                most_common_name = ""
            outfile.write(f"{posid}\t{most_common_name}\n")

    print("集計が完了しました。出力ファイル:", args.output_file)

if __name__ == "__main__":
    main()
