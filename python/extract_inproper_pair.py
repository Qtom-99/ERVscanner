#!/usr/bin/env python3
import argparse

def is_abnormal_pair(rec1, rec2):
    """
    rec1, rec2: BEDファイルの各行をタブ区切りでsplitしたリスト（少なくとも6カラム必要）
      [0] chromosome
      [1] start
      [2] end
      [3] read name
      [4] MAPQ（使用しません）
      [5] strand

    異常ペアとする条件は以下の通り:
      - 同一染色体であること
      - 両リードの開始位置の差が 10,000 塩基以内 (< 10,000)であること
      - 正常なペアの場合、向きは反対で、「+」側のリードが上流（startが小さい）になっているはず。
        これを満たさない場合を異常（abnormal）と判断します。
    """
    # 同一染色体でなければ対象外
    if rec1[0] != rec2[0]:
        return False

    try:
        start1 = int(rec1[1])
        start2 = int(rec2[1])
    except ValueError:
        return False

    # 距離条件: 両リードの開始位置の差が10,000塩基以上なら対象外
    if abs(start1 - start2) >= 10000:
        return False

    strand1 = rec1[5]
    strand2 = rec2[5]

    # 両リードが同じ向きの場合は異常
    if strand1 == strand2:
        return True

    # 向きが反対の場合、正常なペアなら
    # ・もし rec1 が "+"、rec2 が "-" なら、rec1 の start < rec2 の start
    # ・もし rec1 が "-"、rec2 が "+" なら、rec2 の start < rec1 の start
    if strand1 == '+' and strand2 == '-':
        return not (start1 < start2)
    elif strand1 == '-' and strand2 == '+':
        return not (start2 < start1)
    else:
        # 予期しないストランドの値の場合は異常とする
        return True

def extract_abnormal_pairs(input_bed, output_bed):
    with open(input_bed, 'r') as fin, open(output_bed, 'w') as fout:
        # ヘッダー行や空行は除外
        lines = (line for line in fin if not line.startswith("#") and line.strip())
        # 2行ずつ読み込む (/1 のレコードの後に /2 のレコードが来る前提)
        for line1, line2 in zip(lines, lines):
            fields1 = line1.rstrip("\n").split("\t")
            fields2 = line2.rstrip("\n").split("\t")
            if len(fields1) < 6 or len(fields2) < 6:
                continue

            # read name の末尾が /1 と /2 であることを確認し、ベースネームが一致するかチェック
            name1 = fields1[3]
            name2 = fields2[3]
            if not (name1.endswith("/1") and name2.endswith("/2")):
                continue
            if name1[:-2] != name2[:-2]:
                continue

            if is_abnormal_pair(fields1, fields2):
                fout.write("\t".join(fields1) + "\n")
                fout.write("\t".join(fields2) + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="名前順ソート済みBEDから、10kb以内で向きが不自然なペアを抽出します。"
    )
    parser.add_argument("-i", "--input", required=True, help="入力BEDファイルのパス")
    parser.add_argument("-o", "--output", required=True, help="出力BEDファイルのパス")
    args = parser.parse_args()

    extract_abnormal_pairs(args.input, args.output)

if __name__ == "__main__":
    main()
