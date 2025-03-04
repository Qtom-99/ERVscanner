# coding: utf-8

import sys

def compare_and_output(input_file_path, output_file_path):
    with open(input_file_path, 'r', encoding='utf-8') as input_file, \
         open(output_file_path, 'w', encoding='utf-8') as output_file:
        
        lines = input_file.readlines()
        total_lines = len(lines)

        for i in range(0, total_lines, 2):
            if i + 1 < total_lines:
                # 1列目の値が異なる場合には出力ファイルに書き込む
                if lines[i].split('\t')[0] != lines[i + 1].split('\t')[0]:
                    output_file.write(lines[i])
                    output_file.write(lines[i + 1])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    compare_and_output(input_file_path, output_file_path)
