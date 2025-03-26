#!/usr/bin/env python3
import argparse

def is_abnormal_pair(rec1, rec2):
# rec1, rec2: Each is a tab-split list representing a BED file row (requires at least 6 columns)
#   [0] chromosome
#   [1] start
#   [2] end
#   [3] read name
#   [4] MAPQ (not used)
#   [5] strand
#
# Criteria for identifying an abnormal pair:
# - Both reads must be on the same chromosome
# - The difference between their start positions must be less than 10,000 bases
# - In a proper pair, the strands must be opposite, and the "+" strand read should be upstream (i.e., have a smaller start position)
#   If these conditions are not met, the pair is considered abnormal
    
    # Exclude pairs that are not on the same chromosome
    if rec1[0] != rec2[0]:
        return False

    try:
        start1 = int(rec1[1])
        start2 = int(rec2[1])
    except ValueError:
        return False

    # Distance condition: Exclude pairs if the difference between their start positions is 10,000 bases or more
    if abs(start1 - start2) >= 10000:
        return False

    strand1 = rec1[5]
    strand2 = rec2[5]

    # 両リードが同じ向きの場合は異常
    if strand1 == strand2:
        return True

# If the strands are opposite, a proper pair should meet the following conditions:
# - If rec1 is "+" and rec2 is "-", then rec1's start must be less than rec2's start
# - If rec1 is "-" and rec2 is "+", then rec2's start must be less than rec1's start
    if strand1 == '+' and strand2 == '-':
        return not (start1 < start2)
    elif strand1 == '-' and strand2 == '+':
        return not (start2 < start1)
    else:
        # Treat as abnormal if strand values are unexpected
        return True

def extract_abnormal_pairs(input_bed, output_bed):
    with open(input_bed, 'r') as fin, open(output_bed, 'w') as fout:
        # Skip header lines and empty lines
        lines = (line for line in fin if not line.startswith("#") and line.strip())
        # Read two lines at a time (assuming a /1 record is followed by its corresponding /2 record)
        for line1, line2 in zip(lines, lines):
            fields1 = line1.rstrip("\n").split("\t")
            fields2 = line2.rstrip("\n").split("\t")
            if len(fields1) < 6 or len(fields2) < 6:
                continue

            # Ensure that the read names end with /1 and /2, and check that their base names match


            name1 = fields1[3]
            name2 = fields2[3]
            if not (name1.endswith("/1") and name2.endswith("/2")):
                continue
            if name1[:-2] != name2[:-2]:
                continue

            if is_abnormal_pair(fields1, fields2):
                fout.write("\t".join(fields1) + "\n")
                fout.write("\t".join(fields2) + "\n")

def main():
    parser = argparse.ArgumentParser(
        description="# Extract read pairs from a name-sorted BED file that are within 10kb and have abnormal orientation"
    )
    parser.add_argument("-i", "--input", required=True, help="# Path to the input BED file")
    parser.add_argument("-o", "--output", required=True, help="# Path to the input BED file")
    args = parser.parse_args()

    extract_abnormal_pairs(args.input, args.output)

if __name__ == "__main__":
    main()
