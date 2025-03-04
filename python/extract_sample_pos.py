#!/usr/bin/env python3

import sys

def main(sample_name, mutations_tsv, positions_bed, output_bed):
    # サンプルが持つ変異のポジション名を取得
    mutation_positions = set()
    with open(mutations_tsv, 'r') as tsv_file:
        lines = tsv_file.readlines()
        headers = lines[0].strip().split('\t')
        sample_found = False
        for line in lines[1:]:
            fields = line.strip().split('\t')
            if fields[0] == sample_name:
                sample_found = True
                for i in range(1, len(fields)):
                    if fields[i] == '1':
                        mutation_positions.add(headers[i])
                break
        if not sample_found:
            print(f'サンプル "{sample_name}" が {mutations_tsv} に見つかりませんでした。')
            sys.exit(1)

    # BEDファイルから該当するポジションを抽出
    with open(positions_bed, 'r') as bed_file, open(output_bed, 'w') as out_file:
        for line in bed_file:
            fields = line.strip().split('\t')
            if len(fields) >= 4 and fields[3] in mutation_positions:
                out_file.write(line)

    print(f'サンプル "{sample_name}" の変異のみを含むBEDファイル "{output_bed}" を作成しました。')

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('使用法: python script.py サンプル名 mutations.tsv positions.bed 出力ファイル.bed')
        sys.exit(1)
    sample_name = sys.argv[1]
    mutations_tsv = sys.argv[2]
    positions_bed = sys.argv[3]
    output_bed = sys.argv[4]
    main(sample_name, mutations_tsv, positions_bed, output_bed)
