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
cat $DATA_PATH/check_seq/master_file/sampledata/*strand.csv > $DATA_PATH/check_seq/master_file/allsample_data_merge.csv
sort -k5,5 $DATA_PATH/check_seq/master_file/allsample_data_merge.csv > $DATA_PATH/check_seq/master_file/allsample_data_merge_sort.csv
paste $DATA_PATH/check_seq/master_file/allsample_data_merge_sort.csv $DATA_PATH/check_seq/bwa/result/mapping_info_sort > $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master.txt
sort -V -k1,1 -k11,11 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master.txt > $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt
cut -f 1,14 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt > $DATA_PATH/check_seq/master_file/allsample_pos_subject
awk -F'\t' '{split($2, a, "::"); $2 = a[1]; print $1 "\t" $2}' $DATA_PATH/check_seq/master_file/allsample_pos_subject > $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr
python3 $PY_PATH/ERVname_to_class.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr $DATA_PATH/dfam_info/Dfam_ERV_info $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_sed
python3 $PY_PATH/DFaccession_to_class.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_sed $DATA_PATH/dfam_info/Dfam_ERV_info $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_class
python3 $PY_PATH/major_class.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_class $DATA_PATH/check_seq/master_file/pos_allclass_report $DATA_PATH/check_seq/master_file/pos_majorclass_report
awk -F'\t' -v threshold="$IDENTITY_THRESHOLD" '$3 >= threshold' "$DATA_PATH/check_seq/master_file/pos_majorclass_report" > "$DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up"
python3 $PY_PATH/extract_wanted_accession.py $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up $ACCESSIONLIST $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted
cut -f 1 $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted | sort -V > $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted_sort
cut -f 1 $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted | sort -V > $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted_sort
sort -V -k1,1 $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted > $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted_sort_ID
cut -f 1,10 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt > $DATA_PATH/check_seq/master_file/pos_sample
sort -V -k1,1 -k2,2 $DATA_PATH/check_seq/master_file/pos_sample | uniq > $DATA_PATH/check_seq/master_file/pos_sample_uniq
cut -f 1,10,11,13 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt > $DATA_PATH/check_seq/master_file/allsample_POSID_strand_seq
python3 $PY_PATH/get_strand_majoority.py $DATA_PATH/check_seq/master_file/allsample_POSID_strand_seq $DATA_PATH/check_seq/master_file/allsample_POSID_strand_pass $DATA_PATH/check_seq/master_file/allsample_POSID_strand_filtered
python3 $PY_PATH/make01table.py $DATA_PATH/check_seq/master_file/pos_sample_uniq $SAMPLE $DATA_PATH/01table/allsample_01table.tsv
python3 $PY_PATH/extract_major_pos.py $DATA_PATH/01table/allsample_01table.tsv $DATA_PATH/check_seq/master_file/pos_majorclass_${IDENTITY_THRESHOLD//./}up_wanted_sort $DATA_PATH/01table/allsample_01table_${IDENTITY_THRESHOLD//./}up.tsv
date
echo "=== process7: filtering loci by checking the insertion contents and directions: finished extracting loci with ${IDENTITY_THRESHOLD} identity threshold ==="
