import sys

def filter_and_save(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 5列目の値をキーとして行を分類
        rows_by_fifth_column = {}
        for line in infile:
            columns = line.strip().split('\t')
            fifth_column_value = columns[4]

            # 5列目の値をキーとして辞書に行を追加
            if fifth_column_value in rows_by_fifth_column:
                rows_by_fifth_column[fifth_column_value].append(line)
            else:
                rows_by_fifth_column[fifth_column_value] = [line]

        # 5列目の値が同じで、4列目が順に"+"、"-"の行を新しいファイルに書き込む
        for key, rows in rows_by_fifth_column.items():
            # "+"、"-" の順になっているか確認
            plus_minus_values = [row.strip().split('\t')[3] for row in rows]
            if len(plus_minus_values) > 1 and plus_minus_values[0] == "+" and plus_minus_values[1] == "-":
                for row in rows:
                    outfile.write(row)

if __name__ == "__main__":
    # コマンドライン引数の取得
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.tsv output_file.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # 関数を呼び出して行を抽出し、新しいファイルに保存
    filter_and_save(input_file, output_file)

