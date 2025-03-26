import argparse

def create_replacement_dict(dictionary_file):
    replacement_dict = {}
    with open(dictionary_file, 'r', encoding='utf-8') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) >= 3:
                key = columns[1]  # 2nd column (index1)
                value = columns[2]  # 3rd column (index2)
                replacement_dict[key] = value
    return replacement_dict

def replace_second_column_add_third(input_file, replacement_dict, output_file):
    with open(input_file, 'r', encoding='utf-8') as f_in, \
         open(output_file, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            columns = line.strip().split('\t')
            if len(columns) >= 2:  # if second column exists
                original_value = columns[1]
                replaced_value = replacement_dict.get(original_value, original_value)
                columns.append(replaced_value)  # add replaced string to 3rd column
            else:
                columns.append('')  # add space if 2nd column does not exist
            f_out.write('\t'.join(columns) + '\n')  # tab-separated text

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
