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
echo "=== process6: checking ID of insertions ==="
while read line
do
python3 $PY_PATH/find_POS_key4.py $DATA_PATH/sampledata/$line/read_info/${line} $DATA_PATH/allsample_merge/POSITION_LIST $DATA_PATH/sampledata/$line/read_info/${line}_data.csv $NCORE
python3 $PY_PATH/add_strand.py $DATA_PATH/sampledata/$line/read_info/${line}_data.csv $DATA_PATH/sampledata/$line/read_info/${line}_*reads_in_cluster.bed > $DATA_PATH/sampledata/$line/read_info/${line}_data_strand.csv
cp $DATA_PATH/sampledata/$line/read_info/${line}_data_strand.csv $DATA_PATH/check_seq/master_file/sampledata/.
echo "=== process6: checking ID of insertions: completed ${line} ==="
done < $SAMPLE
echo "=== process6: checking ID of insertions: finished ==="
