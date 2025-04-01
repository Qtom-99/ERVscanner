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
            # Check if XA or SA tag exists
            xa_tag = read.get_tag("XA") if read.has_tag("XA") else None
            sa_tag = read.get_tag("SA") if read.has_tag("SA") else None

            # Initialize a flag to indicate if the read should be removed
            remove_read = False

            # Check XA tag if available
            if xa_tag:
                xa_alignments = xa_tag.split(';')[:-1]  # Remove last empty entry
                for aln in xa_alignments:
                    fields = aln.split(',')
                    if len(fields) >= 3:
                        remove_read = True

            # Check SA tag if XA tag didn't trigger removal
            if sa_tag and not remove_read:
                sa_alignments = sa_tag.split(';')[:-1]
                for aln in sa_alignments:
                    fields = aln.split(',')
                    if len(fields) >= 3:
                        chrom = fields[0]
                        cigar = fields[3]  # In SA, CIGAR is at position 4 (0-indexed)
                        if any(chrom == exc for exc in exclude_chroms) and not re.match(r'^\d+S', cigar):
                            remove_read = True
                            break
            # Write to bam_removed if the read meets the criteria, otherwise to bam_out
            if remove_read:
                bam_removed.write(read)
            else:
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
