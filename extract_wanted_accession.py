import sys

def filter_tsv_by_accession(tsv_file, accession_file, output_file):
    # アクセッション番号をセットに読み込む
    with open(accession_file, 'r') as af:
        accession_numbers = set(line.strip() for line in af if line.strip())

    # TSVファイルをフィルタリング
    with open(tsv_file, 'r') as tf, open(output_file, 'w') as of:
        for line in tf:
            fields = line.strip().split("\t")
            if len(fields) >= 2 and fields[1] in accession_numbers:
                of.write(line)

if __name__ == "__main__":
    # 引数の確認
    if len(sys.argv) != 4:
        print("Usage: python filter_tsv.py <tsv_file> <accession_file> <output_file>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    accession_file = sys.argv[2]
    output_file = sys.argv[3]

    filter_tsv_by_accession(tsv_file, accession_file, output_file)
