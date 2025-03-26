#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: {} file1.txt file2.txt\n".format(sys.argv[0]))
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    # # Read file2 and create a dictionary using the 4th column (read name without the part after "/") as the key　and the 6th column (strand) as the value.
    strand_map = {}
    with open(file2_path, "r") as f2:
        for line in f2:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 7:
                continue
            # 例: "ERR019244.4606667/1" → "ERR019244.4606667"
            read_name = cols[3].split('/')[0]
            strand = cols[5]
            strand_map[read_name] = strand

    # Read file1, look up the strand corresponding to the 5th column (key) and append it as the 11th column in the output.
    with open(file1_path, "r") as f1:
        for line in f1:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                print(line)
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                sys.stderr.write("Warning: skipping line with less than 5 columns: {}\n".format(line))
                continue
            # In file1, the key is the 5th column (index 4).
            key = cols[4]
            # If the key exists in file2, set the strand; otherwise, set it to "NA".
            strand = strand_map.get(key, "NA")
            cols.append(strand)
            print("\t".join(cols))

if __name__ == "__main__":
    main()
