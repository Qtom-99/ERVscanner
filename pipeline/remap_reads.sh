#!/bin/bash

# get config file
CONFIG_FILE="$1"

# checking the existence of config file
if [ -z "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not specified."
    echo "Usage: $0 <config_file>"
    exit 1
fi

#  checking the existence of config file
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found."
    exit 1
fi

# reading confing file
eval $(awk -F: '
    /^[[:space:]]*$/ {next}      # skip empty line
    /^[[:space:]]*#/ {next}      # ignore comment line

    # remove spaces
    {gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)}

    # exporting parameters
    {print $1"=\"" $2 "\""}
' "$CONFIG_FILE")

TARGET_REPEAT_FASTA="$DATA_PATH/check_seq/bwa/subject/nrph.fasta"

date

echo "=== process4: merging fastq ==="
cat $DATA_PATH/check_seq/fq/sampledata/*.fq > $DATA_PATH/check_seq/fq/allsample_ERVread.fq

echo "=== remap to Dfam (keep secondary + supplementary) ==="
bwa mem -t "$NCORE" -a "$TARGET_REPEAT_FASTA" \
  "$DATA_PATH/check_seq/fq/allsample_ERVread.fq" \
| samtools view -@ "$NCORE" -b -o \
  "$DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV.all.bam"

samtools view -@ $NCORE -F 2048 -o $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs_F2048.bam $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV.all.bam

samtools view $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs_F2048.bam | cut -f 1-6 - | sort -k1,1 - > $DATA_PATH/check_seq/bwa/result/mapping_info_sort

samtools sort -@ "$NCORE" -n \
  -o "$DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV.all.qnamesort.bam" \
  "$DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV.all.bam"

samtools view "$DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV.all.qnamesort.bam" \
| awk 'BEGIN{OFS="\t"}
{
  q=$1; flag=$2; r=$3; pos=$4; mapq=$5; cigar=$6;

  # secondary/supplementary
  is_sec = and(flag,256)?1:0;
  is_supp = and(flag,2048)?1:0;

  AS="."; NM=".";
  for(i=12;i<=NF;i++){
    if($i ~ /^AS:i:/){AS=$i; sub(/^AS:i:/,"",AS)}
    else if($i ~ /^NM:i:/){NM=$i; sub(/^NM:i:/,"",NM)}
  }

  print q, flag, is_sec, is_supp, r, pos, mapq, cigar, AS, NM
}' \
> "$DATA_PATH/check_seq/bwa/result/mapping_info_all"

echo "=== process4: finishing fastq merge and remapping ==="
date
echo "=== process5: making a list of insertion positions ==="
cat $DATA_PATH/allsample_merge/sampledata/*.bed > $DATA_PATH/allsample_merge/allsample_merged.bed
sort -V -k1,1 -k2,2 $DATA_PATH/allsample_merge/allsample_merged.bed | bedtools merge -i - -d 50 > $DATA_PATH/allsample_merge/allsample_merged_d50.bed
awk 'BEGIN {OFS="\t"} {print "POS" NR, $1, $2, $3}' $DATA_PATH/allsample_merge/allsample_merged_d50.bed > $DATA_PATH/allsample_merge/POSITION_LIST
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "POS" NR}' $DATA_PATH/allsample_merge/allsample_merged_d50.bed > $DATA_PATH/allsample_merge/pos.bed
date
echo "=== process5: finishing to make a list of insertion positions ==="
