import argparse

def create_replacement_dict(dictionary_file):
    replacement_dict = {}
    with open(dictionary_file, 'r', encoding='utf-8') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) >= 3:
                key = columns[1]  # 2列目（インデックス1）
                value = columns[2]  # 3列目（インデックス2）
                replacement_dict[key] = value
    return replacement_dict

def replace_second_column_add_third(input_file, replacement_dict, output_file):
    with open(input_file, 'r', encoding='utf-8') as f_in, \
         open(output_file, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            columns = line.strip().split('\t')
            if len(columns) >= 2:  # 2列目が存在する場合
                original_value = columns[1]
                replaced_value = replacement_dict.get(original_value, original_value)
                columns.append(replaced_value)  # 置換後の値を3列目として追加
            else:
                columns.append('')  # 2列目がない場合は空白を追加
            f_out.write('\t'.join(columns) + '\n')  # タブ区切りで書き込み

def main():
    parser = argparse.ArgumentParser(description='Replace the second column in a text file and append the result as the third column.')
    parser.add_argument('input_file', help='Path to the input text file.')
    parser.add_argument('dictionary_file', help='Path to the dictionary text file.')
    parser.add_argument('output_file', help='Path to the output text file.')

    args = parser.parse_args()

    replacement_dict = create_replacement_dict(args.dictionary_file)
    replace_second_column_add_third(args.input_file, replacement_dict, args.output_file)

if __name__ == '__main__':
    main()
