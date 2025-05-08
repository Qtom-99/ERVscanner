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

date
echo "=== process9: storing the result of genotyping  ==="
python3 $PY_PATH/renew_genotype_0123.py $DATA_PATH/01table/allsample_01table_07up.tsv $DATA_PATH/genotype/allsample $DATA_PATH/01table/genotype/allsample_07_sort_genotype.tsv
python3 $PY_PATH/remove_all0pos_strand.py $DATA_PATH/01table/genotype/allsample_07_sort_genotype.tsv $DATA_PATH/check_seq/master_file/allsample_POSID_strand_filtered $DATA_PATH/01table/genotype/allsample_07_sort_genotype_remove0.tsv
date
echo "=== process10: estimating breakpoints  ==="
cat $DATA_PATH/genotype/allsample/*.txt | sort -V -k4,4 - > $DATA_PATH/genotype/allsample_genotype_sort.txt
python3 $PY_PATH/get_pos_breakpoint.py $DATA_PATH/genotype/allsample_genotype_sort.txt $DATA_PATH/genotype/allsample_breakpoint.txt
sort -V -k1,1 -k2,2 $DATA_PATH/genotype/allsample_breakpoint.txt > $DATA_PATH/genotype/allsample_breakpoint_sort.txt
date
echo "=== process11: making vcf file ==="
python3 $PY_PATH/get_ervname_majority.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_class $DATA_PATH/check_seq/master_file/allsample_pos_ervname_majority
python3 $PY_PATH/make_vcf_multcore.py --file1 $DATA_PATH/genotype/allsample_breakpoint_sort.txt --file2 $DATA_PATH/allsample_merge/pos.bed --file3 $DATA_PATH/01table/genotype/allsample_07_sort_genotype_remove0.tsv --file4 $SAMPLE --file5 $DATA_PATH/check_seq/master_file/pos_majorclass_07up_wanted_sort_ID --file6 $DATA_PATH/check_seq/master_file/allsample_pos_ervname_majority --file7 $DATA_PATH/check_seq/master_file/allsample_POSID_strand_pass --file8 $DATA_PATH/genotype/allsample --output $DATA_PATH/vcf/allsample_ervins.vcf
