#!/bin/bash

# default

# getting option values
while getopts "d:s:" opt; do
  case $opt in
    s) SAMPLE="$OPTARG" ;;
    d) DATA_PATH="$OPTARG" ;;
    \?) echo "Usage: $0 [-i input] [-o output]" >&2; exit 1 ;;
  esac
done

if [[ -z "$SAMPLE" || -z "$DATA_PATH" ]]; then
  echo "Error: requied option values are missing" >&2
  exit 1
fi

date
echo "=== process9: storing the result of genotyping  ==="
python3 $DATA_PATH/script/renew_genotype_0123.py $DATA_PATH/01table/allsample_01table_07up.tsv $DATA_PATH/genotype/allsample $DATA_PATH/01table/genotype/allsample_07_sort_genotype.tsv
python3 $DATA_PATH/script/remove_all0pos_strand.py $DATA_PATH/01table/genotype/allsample_07_sort_genotype.tsv $DATA_PATH/check_seq/master_file/allsample_POSID_strand_filtered $DATA_PATH/01table/genotype/allsample_07_sort_genotype_remove0.tsv
date
echo "=== process10: estimating breakpoints  ==="
cat $DATA_PATH/genotype/allsample/*.txt | sort -V -k4,4 - > $DATA_PATH/genotype/allsample_genotype_sort.txt
python3 $DATA_PATH/script/get_pos_breakpoint.py $DATA_PATH/genotype/allsample_genotype_sort.txt $DATA_PATH/genotype/allsample_breakpoint.txt
sort -V -k1,1 -k2,2 $DATA_PATH/genotype/allsample_breakpoint.txt > $DATA_PATH/genotype/allsample_breakpoint_sort.txt
date
echo "=== process11: making vcf file ==="
python3 $DATA_PATH/script/get_ervname_majority.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_class $DATA_PATH/check_seq/master_file/allsample_pos_ervname_majority
python3 $DATA_PATH/script/make_vcf_multcore.py --file1 $DATA_PATH/genotype/allsample_breakpoint_sort.txt --file2 $DATA_PATH/allsample_merge/pos.bed --file3 $DATA_PATH/01table/genotype/allsample_07_sort_genotype_remove0.tsv --file4 $SAMPLE --file5 $DATA_PATH/check_seq/master_file/pos_majorclass_07up_wanted_sort_ID --file6 $DATA_PATH/check_seq/master_file/allsample_pos_ervname_majority --file7 $DATA_PATH/check_seq/master_file/allsample_POSID_strand_pass --file8 $DATA_PATH/genotype/allsample --output $DATA_PATH/vcf/allsample_ervins.vcf
