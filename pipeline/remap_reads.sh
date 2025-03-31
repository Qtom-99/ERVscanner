#!/bin/bash

# default values
NCORE=1

# getting option values
while getopts "d:n:" opt; do
  case $opt in
    d) DATA_PATH="$OPTARG" ;;
    \?) echo "Usage: $0 [-i input] [-o output]" >&2; exit 1 ;;
  esac
done

TARGET_REPEAT_FASTA="$DATA_PATH/check_seq/bwa/subject/nrph.fasta"

if [[ -z "$DATA_PATH" ]]; then
  echo "Error: requied option values are missing" >&2
  exit 1
fi

if [[ ! -d "$DATA_PATH" ]]; then
  echo "Error: Directory '$DATA_PATH' not found." >&2
  exit 1
fi

date

echo "=== process4: merging fastq ==="
cat $DATA_PATH/check_seq/fq/sampledata/*.fq > $DATA_PATH/check_seq/fq/allsample_ERVread.fq
bwa mem -t $NCORE -o $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs.bam $TARGET_REPEAT_FASTA $DATA_PATH/check_seq/fq/allsample_ERVread.fq
samtools view -@ $NCORE -F 2048 -o $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs_F2048.bam $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs.bam
samtools view $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs_F2048.bam | cut -f 1-6 - | sort -k1,1 - > $DATA_PATH/check_seq/bwa/result/mapping_info_sort
echo "=== process4: finishing fastq merge and remapping ==="
date
echo "=== process5: making a list of insertion positions ==="
cat $DATA_PATH/allsample_merge/sampledata/*.bed > $DATA_PATH/allsample_merge/allsample_merged.bed
sort -V -k1,1 -k2,2 $DATA_PATH/allsample_merge/allsample_merged.bed | bedtools merge -i - -d 50 > $DATA_PATH/allsample_merge/allsample_merged_d50.bed
awk 'BEGIN {OFS="\t"} {print "POS" NR, $1, $2, $3}' $DATA_PATH/allsample_merge/allsample_merged_d50.bed > $DATA_PATH/allsample_merge/POSITION_LIST
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "POS" NR}' $DATA_PATH/allsample_merge/allsample_merged_d50.bed > $DATA_PATH/allsample_merge/pos.bed
date
echo "=== process5: finishing to make a list of insertion positions ==="
