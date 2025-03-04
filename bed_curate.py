import sys

def read_bed_file(file_path):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            read_name = fields[3].split('/')[0]
            data[read_name] = line
    return data

def merge_bed_files(file1_path, file2_path, output_path):
    data1 = read_bed_file(file1_path)
    data2 = read_bed_file(file2_path)

    with open(output_path, 'w') as output_file:
        for read_name in set(data1.keys()) & set(data2.keys()):
            output_file.write(data1[read_name])
            output_file.write(data2[read_name])

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py bed_file1.bed bed_file2.bed bed_file3.bed")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_path = sys.argv[3]

    merge_bed_files(file1_path, file2_path, output_path)
