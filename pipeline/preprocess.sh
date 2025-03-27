#!/bin/bash

# default values

# getting option values
while getopts "i:s:r:t:d:n:b:q:c:a:p:" opt; do
  case $opt in
    s) SAMPLE="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    d) DATA_PATH="$OPTARG" ;;
    f) DFAM_INFO="$OPTARG" ;;
    n) NRPH="$OPTARG" ;;
    t) TARGET="$OPTARG" ;;
    a) ALT_CHR_LIST="$OPTARG" ;;
    p) PY_PATH="$OPTARG" ;;
    \?) echo "Usage: $0 [-i input] [-o output]" >&2; exit 1 ;;
  esac
done

# checking all required options are specified
if [[ -z "$SAMPLE" || -z "$REF_GENOME" || -z "$DATA_PATH" || -z "$DFAM_INFO" || -z "$ALT_CHR_LIST" || -z "$NRPH" || -z "$TARGET" || -z "$PY_PATH" ]]; then
  echo "Error: requied option values are missing" >&2
  exit 1
fi

# chekcing file does exist
if [[ ! -f "$SAMPLE" ]]; then
  echo "Error: File '$SAMPLE' not found." >&2
  exit 1
fi

if [[ ! -f "$REF_GENOME" ]]; then
  echo "Error: File '$REF_GENOME' not found." >&2
  exit 1
fi

if [[ ! -d "$DATA_PATH" ]]; then
  echo "Error: Directory '$DATA_PATH' not found." >&2
  exit 1
fi

if [[ ! -f "$NRPH" ]]; then
  echo "Error: File '$NRPH' not found." >&2
  exit 1
fi

if [[ ! -f "$NRPH" ]]; then
  echo "Error: File '$NRPH' not found." >&2
  exit 1
fi

if [[ ! -f "$TARGET" ]]; then
  echo "Error: File '$TARGET' not found." >&2
  exit 1
fi

if [[ ! -f "$DFAM_INFO" ]]; then
  echo "Error: File '$DFAM_INFO' not found." >&2
  exit 1
fi
if [[ ! -f "$ALT_CHR_LIST" ]]; then
  echo "Error: File '$ALT_CHR_LIST' not found." >&2
  exit 1
fi

if [[ ! -d "$PY_PATH" ]]; then
  echo "Error: Directory '$PY_PATH' not found." >&2
  exit 1
fi

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
mkdir $DATA_PATH/inMR
mkdir $DATA_PATH/inMR/sampledata
mkdir $DATA_PATH/inMR/script
mkdir $DATA_PATH/inMR/check_seq
mkdir $DATA_PATH/inMR/check_seq/fq
mkdir $DATA_PATH/inMR/check_seq/fq/sampledata
mkdir $DATA_PATH/inMR/check_seq/bwa
mkdir $DATA_PATH/inMR/check_seq/bwa/result
mkdir $DATA_PATH/inMR/check_seq/bwa/subject
mkdir $DATA_PATH/inMR/check_seq/master_file
mkdir $DATA_PATH/inMR/check_seq/master_file/sampledata
mkdir $DATA_PATH/inMR/dfam_info
mkdir $DATA_PATH/inMR/allsample_merge
mkdir $DATA_PATH/inMR/allsample_merge/sampledata
mkdir $DATA_PATH/inMR/genotype
mkdir $DATA_PATH/inMR/genotype/breakpoint
mkdir $DATA_PATH/inMR/genotype/breakpoint/sampledata
mkdir $DATA_PATH/inMR/genotype/sampledata
mkdir $DATA_PATH/inMR/genotype/allsample
mkdir $DATA_PATH/inMR/01table
mkdir $DATA_PATH/inMR/01table/genotype
mkdir $DATA_PATH/inMR/vcf
mkdir $DATA_PATH/reference
while read line
do
mkdir $DATA_PATH/sampledata/$line
mkdir $DATA_PATH/sampledata/$line/clip
mkdir $DATA_PATH/sampledata/$line/read_info
mkdir $DATA_PATH/sampledata/$line/curated_bam
mkdir $DATA_PATH/genotype/sampledata/$line
mkdir $DATA_PATH/inMR/sampledata/$line
mkdir $DATA_PATH/inMR/sampledata/$line/clip
mkdir $DATA_PATH/inMR/sampledata/$line/read_info
mkdir $DATA_PATH/inMR/sampledata/$line/curated_bam
mkdir $DATA_PATH/inMR/genotype/sampledata/$line
done < $SAMPLE
echo "=== copying and generating files ==="
cp $PY_PATH/*.py $DATA_PATH/script/ 
cp $REF_GENOME $DATA_PATH/reference/reference.fasta
cp $DFAM_INFO $DATA_PATH/dfam_info/
cut -f2 $DFAM_INFO> | uniq | sort > $DATA_PTH/dfam_info/target_name.txt
cut -f3 $DFAM_INFO> | uniq | sort > $DATA_PTH/dfam_info/target_class.txt
python3 $DATA_PATH/script/wordgrep.py $TARGET $DATA_PTH/dfam_info/target_name.txt $DATA_PATH/dfam_info/target.hits.gz
python3 $DATA_PATH/script/make_bed.py $DATA_PATH/dfam_info/target.hits.gz $DATA_PATH/dfam_info/target.bed
python3 $DATA_PATH/script/make_bed.py $NRPH $DATA_PATH/dfam_info/nrph.bed
sort -V -k1,1 -k2,2 $DATA_PATH/dfam_info/nrph.bed | bedtools merge -d 50 - > $DATA_PATH/dfam_info/nrph.sorted.bed
bedtools getfasta -fi $DATA_PATH/reference/reference.fasta -name+ -bed $DATA_PATH/dfam_info/nrph.sorted.bed -s -fo $DATA_PATH/check_seq/bwa/subject/nrph.fasta
bwa index $DATA_PATH/check_seq/bwa/subject/nrph.fasta
