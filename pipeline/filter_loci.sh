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

ACCESSIONLIST="$DATA_PATH/dfam_info/target_class.txt"

date
echo "=== process7: filtering loci by checking the insertion contents and directions  ==="
cat $DATA_PATH/check_seq/master_file/sampledata/*strand.tsv > $DATA_PATH/check_seq/master_file/allsample_data_merge.tsv
sort -k5,5 $DATA_PATH/check_seq/master_file/allsample_data_merge.csv > $DATA_PATH/check_seq/master_file/allsample_data_merge_sort.csv
paste $DATA_PATH/check_seq/master_file/allsample_data_merge_sort.csv $DATA_PATH/check_seq/bwa/result/mapping_info_sort > $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master.txt
sort -V -k1,1 -k11,11 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master.txt > $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt
python3 $PY_PATH/identify_family.py $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV.all.qnamesort.bam $DATA_PATH/check_seq/master_file/allsample_data_merge.tsv $DATA_PATH/check_seq/master_file/repeat_per_POS_info.tsv
sort -V -k1,1 $DATA_PATH/check_seq/master_file/repeat_per_POS_info.call_pos.txt > $DATA_PATH/check_seq/master_file/repeat_per_POS_info.call_pos_sorted.txt
awk -F'\t' ' NR==1 {next} $2!="NA" { split($2, a, "::"); print $1 "\t" a[1] } ' $DATA_PATH/check_seq/master_file/repeat_per_POS_info.call.tsv > $DATA_PATH/check_seq/master_file/allsample_ervname.tsv
awk -F'\t' '{print $1 "\t" $5}' $DATA_PATH/check_seq/master_file/repeat_per_POS_info.call.tsv > $DATA_PATH/check_seq/master_file/allsample_strand.tsv
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{map[$2]=$3; next} {if($2 in map){$2=map[$2]} else {$2="NA"}; print}' "$DFAM_INFO" "$DATA_PATH/check_seq/master_file/allsample_ervname.tsv" > "$DATA_PATH/check_seq/master_file/allsample_class.tsv"
cut -f 1,10 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt > $DATA_PATH/check_seq/master_file/pos_sample
sort -V -k1,1 -k2,2 $DATA_PATH/check_seq/master_file/pos_sample | uniq > $DATA_PATH/check_seq/master_file/pos_sample_uniq
awk 'NR==FNR {a[$1]; next} ($1 in a)' $DATA_PATH/check_seq/master_file/repeat_per_POS_info.call_pos_sorted.txt $DATA_PATH/check_seq/master_file/pos_sample_uniq > $DATA_PATH/check_seq/master_file/pos_sample_uniq_filtered
python3 $PY_PATH/make01table.py $DATA_PATH/check_seq/master_file/pos_sample_uniq $SAMPLE $DATA_PATH/01table/allsample_01table.tsv
python3 $PY_PATH/extract_major_pos.py $DATA_PATH/01table/allsample_01table.tsv $DATA_PATH/check_seq/master_file/repeat_per_POS_info.call_pos_sorted.txt $DATA_PATH/01table/allsample_01table_filtered.tsv
date
echo "=== process7: filtering loci by checking the insertion contents and directions: finished extracting loci ==="
