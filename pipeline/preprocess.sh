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

echo "=== creating directories ==="
mkdir $DATA_PATH/sampledata
mkdir $DATA_PATH/script
mkdir $DATA_PATH/check_seq
mkdir $DATA_PATH/check_seq/fq
mkdir $DATA_PATH/check_seq/fq/sampledata
mkdir $DATA_PATH/check_seq/bwa
mkdir $DATA_PATH/check_seq/bwa/result
mkdir $DATA_PATH/check_seq/bwa/subject
mkdir $DATA_PATH/check_seq/master_file
mkdir $DATA_PATH/check_seq/master_file/sampledata
mkdir $DATA_PATH/dfam_info
mkdir $DATA_PATH/allsample_merge
mkdir $DATA_PATH/allsample_merge/sampledata
mkdir $DATA_PATH/genotype
mkdir $DATA_PATH/genotype/breakpoint
mkdir $DATA_PATH/genotype/breakpoint/sampledata
mkdir $DATA_PATH/genotype/sampledata
mkdir $DATA_PATH/genotype/allsample
mkdir $DATA_PATH/01table
mkdir $DATA_PATH/01table/genotype
mkdir $DATA_PATH/vcf
mkdir $DATA_PATH/reference
while read line
do
mkdir $DATA_PATH/sampledata/$line
mkdir $DATA_PATH/sampledata/$line/clip
mkdir $DATA_PATH/sampledata/$line/read_info
mkdir $DATA_PATH/sampledata/$line/curated_bam
mkdir $DATA_PATH/genotype/sampledata/$line
done < $SAMPLE

echo "=== copying and generating files ==="
cp $REF_GENOME $DATA_PATH/reference/reference.fasta
samtools faidx $DATA_PATH/reference/reference.fasta
#cp $DFAM_INFO $DATA_PATH/dfam_info/
#cp $ALT_CHR_LIST $DATA_PATH/reference/alt_chr_list
cut -f2 $DFAM_INFO | uniq | sort > $DATA_PATH/dfam_info/target_name.txt
cut -f3 $DFAM_INFO | uniq | sort > $DATA_PATH/dfam_info/target_class.txt
echo "=== generating target bed file ==="
python3 $PY_PATH/wordgrep.py $ALL_REPEAT $DATA_PATH/dfam_info/target_name.txt $DATA_PATH/dfam_info/target.hits.gz
echo "=== generating target bed file: finished ==="
zcat $DATA_PATH/dfam_info/target.hits.gz | cut -f 1,10,11  > $DATA_PATH/dfam_info/target.hits_convert.bed
python3 $PY_PATH/convert_bed.py $DATA_PATH/dfam_info/target.hits_convert.bed  $DATA_PATH/dfam_info/target_pre.bed 0
sort -V -k1,1 -k2,2 $DATA_PATH/dfam_info/target_pre.bed | bedtools merge -d 50 -i - > $DATA_PATH/dfam_info/target.bed
rm $DATA_PATH/dfam_info/target_pre.bed
zcat $NRPH_REPEAT | awk 'BEGIN {OFS="\t"} {print $1, $10, $11, $3, $3, $9}'  > $DATA_PATH/dfam_info/nrph.hits_convert.bed
python3 $PY_PATH/convert_bed.py  $DATA_PATH/dfam_info/nrph.hits_convert.bed  $DATA_PATH/dfam_info/nrph.bed 1
sort -V -k1,1 -k2,2 $DATA_PATH/dfam_info/nrph.bed  > $DATA_PATH/dfam_info/nrph.sorted.bed
bedtools getfasta -fi $DATA_PATH/reference/reference.fasta -name+ -bed $DATA_PATH/dfam_info/nrph.sorted.bed -s -fo $DATA_PATH/check_seq/bwa/subject/nrph.fasta
bwa index $DATA_PATH/check_seq/bwa/subject/nrph.fasta

