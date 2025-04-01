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

REF_GENOME="$DATA_PATH/reference/reference.fasta"

while read line
do
date
echo "=== process8: SAMPLE:${line} genotyping  ==="
samtools view -h $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.bam | awk '/^@/ || $6 ~ /[0-9]+[SH]/' | samtools view -b -o $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads.bam -
bedtools bamtobed -i $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads.bam > $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads.bed
sort -V -k1,1 -k2,2 $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads.bed > $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort.bed
bedtools merge -i $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort.bed -d 50 -s -c 6 -o distinct > $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50.bed
bedtools cluster -d 50 -i $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50.bed > $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50.bed
python3 $PY_PATH/bed_cluster_plusminus.py $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50.bed $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50_plusminus.bed
python3 $PY_PATH/bed_cluster_plusminus_jun.py $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50_plusminus.bed $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50_plusminus_jun.bed
python3 $PY_PATH/get_center_pos.py $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50_plusminus_jun.bed $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50_plusminus_jun_centerpos.txt
python3 $PY_PATH/extract_sample_pos.py ${line} $DATA_PATH/01table/allsample_01table_${IDENTITY_THRESHOLD//./}up.tsv $DATA_PATH/allsample_merge/pos.bed $DATA_PATH/sampledata/${line}/clip/${line}_pos.bed
python3 $PY_PATH/match_centerpos_to_POS.py $DATA_PATH/sampledata/${line}/clip/${line}_clipped_reads_sort_merged50_clustered50_plusminus_jun_centerpos.txt $DATA_PATH/sampledata/${line}/clip/${line}_pos.bed $DATA_PATH/sampledata/${line}/clip/${line}_breakpoint.tsv
samtools view -@ "$NCORE" -T "$REF_GENOME" -L $DATA_PATH/sampledata/${line}/clip/${line}_pos.bed -o "$DATA_PATH"/genotype/sampledata/${line}/${line}_muttype.bam -O BAM "$INPUT_PATH"/${line}."$INPUT_TYPE"
samtools view -@ $NCORE -q ${MAPQ_THRESHOLD} -o "$DATA_PATH"/genotype/sampledata/${line}/${line}_muttype_q${MAPQ_THRESHOLD}.bam "$DATA_PATH"/genotype/sampledata/${line}/${line}_muttype.bam
samtools view -@ $NCORE -h $DATA_PATH/genotype/sampledata/${line}/${line}_muttype_q${MAPQ_THRESHOLD}.bam | awk '$6 !~ /[SH]/ || $1 ~ /^@/' | samtools view -b -o $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping.bam -
bedtools bamtobed -i $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping.bam > $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping.bed
sort -V -k1,1 -k2,2 $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping.bed > $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping_sort.bed
cut -f 1,2,3 $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping_sort.bed > $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping_sort123.bed
awk -v OFS="\t" '{print $1, $2-15, $2+15, $4}' $DATA_PATH/sampledata/${line}/clip/${line}_breakpoint.tsv > $DATA_PATH/sampledata/${line}/clip/${line}_breakpoint.bed
bedtools intersect -a $DATA_PATH/sampledata/${line}/clip/${line}_breakpoint.bed -b $DATA_PATH/genotype/sampledata/${line}/${line}_q${MAPQ_THRESHOLD}_no_clipping_sort123.bed -f 1.0 -c > $DATA_PATH/genotype/sampledata/$line/${line}_n_of_covreads_noclip.txt
cp $DATA_PATH/genotype/sampledata/$line/${line}_n_of_covreads_noclip.txt $DATA_PATH/genotype/allsample/.
date
echo "=== process8: SAMPLE:${line} genotyping: finished  ==="
done < $SAMPLE
