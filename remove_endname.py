# coding: utf-8

import sys

def remove_suffix(input_file_path, output_file_path):
    with open(input_file_path, 'r') as input_file, \
         open(output_file_path, 'w') as output_file:
        
        for line in input_file:
            # 行から"/1"および"/2"を削除して出力ファイルに書き込む
            modified_line = line.replace("/1", "").replace("/2", "")
            output_file.write(modified_line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    remove_suffix(input_file_path, output_file_path)
