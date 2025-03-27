#!/bin/bash

# default values
NCORE=1

# getting option values
while getopts "d:s:n:" opt; do
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

# chekcing file does exist
if [[ ! -f "$SAMPLE" ]]; then
  echo "Error: File '$SAMPLE' not found." >&2
  exit 1
fi

if [[ ! -d "$DATA_PATH" ]]; then
  echo "Error: Directory '$DATA_PATH' not found." >&2
  exit 1
fi

date
echo "=== process6: checking ID of insertions ==="
while read line
do
python3 $DATA_PATH/script/find_POS_key4.py $DATA_PATH/sampledata/$line/read_info/${line} $DATA_PATH/allsample_merge/POSITION_LIST $DATA_PATH/sampledata/$line/read_info/${line}_data.csv $NCORE
python3 $DATA_PATH/script/add_strand.py $DATA_PATH/sampledata/$line/read_info/${line}_data.csv $DATA_PATH/sampledata/$line/read_info/${line}_*reads_in_cluster.bed > $DATA_PATH/sampledata/$line/read_info/${line}_data_strand.csv
cp $DATA_PATH/sampledata/$line/read_info/${line}_data_strand.csv $DATA_PATH/check_seq/master_file/sampledata/.
echo "=== process6: checking ID of insertions: completed ${line} ==="
done < $SAMPLE
echo "=== process6: checking ID of insertions: finished ==="
