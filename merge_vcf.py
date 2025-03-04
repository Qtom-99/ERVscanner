#!/usr/bin/env python3
import argparse
import re

def read_header(vcf_file):
    """ヘッダー行（##～および#CHROM行）を返す"""
    header_lines = []
    with open(vcf_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
            else:
                # ヘッダーの終わりに達したらbreak
                break
    return header_lines

def process_file1(vcf_file):
    """ファイル1のレコードを読み込み、INFOにMI=-を追加"""
    records = []
    with open(vcf_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            # INFOフィールドは8列目（0-indexedで7番目）
            if fields[7] == "." or fields[7] == "":
                fields[7] = "MI=-"
            else:
                fields[7] = fields[7] + ";MI=-"
            records.append(fields)
    return records

def process_file2(vcf_file):
    """
    ファイル2のレコードを読み込み、INFOにMI=+を追加し、
    ID（3列目）の数字部分を保持して、MR<数字>に置換する。
    例: 元のIDが"POS25"なら"MR25"とする。
    """
    records = []
    with open(vcf_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[7] == "." or fields[7] == "":
                fields[7] = "MI=+"
            else:
                fields[7] = fields[7] + ";MI=+"
            original_id = fields[2]
            m = re.search(r'(\d+)', original_id)
            if m:
                # 数字部分を保持してMR<数字>に置換
                fields[2] = f"MR{m.group(1)}"
            else:
                # 数字が抽出できない場合は、元のIDにMRを付与
                fields[2] = f"MR{original_id}"
            records.append(fields)
    return records

def sort_records(records):
    """CHROMとPOSをキーにしてレコードをソート（POSは数値として扱う）"""
    def sort_key(fields):
        chrom = fields[0]
        try:
            pos = int(fields[1])
        except ValueError:
            pos = 0
        return (chrom, pos)
    return sorted(records, key=sort_key)

def add_MI_header(header_lines):
    """ヘッダーにMIの定義が無ければ追加する（#CHROM行の直前に挿入）"""
    mi_defined = any("##INFO=<ID=MI," in line for line in header_lines)
    if mi_defined:
        return header_lines
    new_header = []
    for line in header_lines:
        if line.startswith("#CHROM"):
            new_header.append('##INFO=<ID=MI,Number=1,Type=String,Description="masked region insertion">')
        new_header.append(line)
    return new_header

def write_merged_vcf(header_lines, records, output_file):
    with open(output_file, "w") as out:
        for line in header_lines:
            out.write(line + "\n")
        for fields in records:
            out.write("\t".join(fields) + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="2つのVCFファイルをマージします。ファイル1はINFOにMI=-、"
                    "ファイル2はINFOにMI=+を追加し、IDは元の数字部分を保持してMR<ID>に置換します。"
    )
    parser.add_argument("-i1", "--input1", required=True, help="入力VCFファイル1のパス")
    parser.add_argument("-i2", "--input2", required=True, help="入力VCFファイル2のパス")
    parser.add_argument("-o", "--output", required=True, help="出力マージVCFファイルのパス")
    args = parser.parse_args()

    file1 = args.input1
    file2 = args.input2
    output_file = args.output

    # ヘッダーはファイル1から読み込み、MIタグ定義を追加
    header_lines = read_header(file1)
    header_lines = add_MI_header(header_lines)

    # レコードの処理
    records1 = process_file1(file1)
    records2 = process_file2(file2)

    # 両ファイルのレコードを結合してソート（CHROM, POS順）
    all_records = records1 + records2
    sorted_records = sort_records(all_records)

    # マージVCFの書き出し
    write_merged_vcf(header_lines, sorted_records, output_file)
    print(f"マージ完了: {output_file}")

if __name__ == "__main__":
    main()
