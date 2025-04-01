import pysam
import argparse
import re

def load_exclude_chroms(file_path):
    """
    Load a list of chromosomes to exclude from a text file (one per line).
    :param file_path: Path to the text file containing chromosome names.
    :return: List of chromosome names to exclude.
    """
    with open(file_path, "r") as f:
        return [line.strip() for line in f if line.strip()]

def filter_bam_by_xa(input_bam, output_bam, removed_bam, exclude_chroms):
    """
    Filter reads based on the XA tag in a BAM file.
    Reads containing any excluded chromosome in the XA tag will be removed 
    and saved in a separate BAM file.

    :param input_bam: Path to the input BAM file.
    :param output_bam: Path to the filtered BAM file (retained reads).
    :param removed_bam: Path to the BAM file storing removed reads.
    :param exclude_chroms: List of chromosome names to filter out.
    """
    with pysam.AlignmentFile(input_bam, "rb") as bam_in, \
         pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out, \
         pysam.AlignmentFile(removed_bam, "wb", header=bam_in.header) as bam_removed:
        
        for read in bam_in:
            xa_tag = read.get_tag("XA") if read.has_tag("XA") or read.has_tag("SA") else None
            cigar = xa_tag.split(':')[2].split(',')[4]
            if xa_tag and any(chrom in xa_tag for chrom in exclude_chroms) and re.match('r/\d+M/', cigar):
                # Write removed reads to a separate BAM file
                bam_removed.write(read)
            else:
                # Write retained reads to the output BAM file
                bam_out.write(read)

def main():
    parser = argparse.ArgumentParser(description="Filter reads by XA tag and save removed reads separately.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file (filtered reads)")
    parser.add_argument("removed_bam", help="BAM file to store removed reads")
    parser.add_argument("exclude_list", help="File containing a list of chromosomes to exclude (one per line)")

    args = parser.parse_args()
    
    exclude_chroms = load_exclude_chroms(args.exclude_list)
    
    filter_bam_by_xa(args.input_bam, args.output_bam, args.removed_bam, exclude_chroms)

if __name__ == "__main__":
    main()
