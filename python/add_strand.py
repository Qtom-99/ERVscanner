#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: {} file1.txt file2.txt\n".format(sys.argv[0]))
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    # ファイル２を読み込み、4列目のリード名（"/"以降を除去）をキーに、6列目（strand）を辞書に登録
    strand_map = {}
    with open(file2_path, "r") as f2:
        for line in f2:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 7:
                continue
            # 例: "ERR019244.4606667/1" → "ERR019244.4606667"
            read_name = cols[3].split('/')[0]
            strand = cols[5]
            strand_map[read_name] = strand

    # ファイル１を読み込み、5列目（キー）に対応するstrandを11列目として追加して出力
    with open(file1_path, "r") as f1:
        for line in f1:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                print(line)
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                sys.stderr.write("Warning: skipping line with less than 5 columns: {}\n".format(line))
                continue
            # ファイル１では5列目（インデックス4）がキー
            key = cols[4]
            # ファイル２にキーがあればstrand、なければ "NA" を設定
            strand = strand_map.get(key, "NA")
            cols.append(strand)
            print("\t".join(cols))

if __name__ == "__main__":
    main()
