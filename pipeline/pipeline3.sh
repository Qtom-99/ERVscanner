#!/bin/bash

SAMPLE=$1
DATA_PATH=$2
NCORE=$3

date
echo "=== process6: 各サンプルのinsertionとポジションIDを照合中 ==="
while read line
do
python3 $DATA_PATH/script/find_POS_key4.py $DATA_PATH/sampledata/$line/read_info/${line} $DATA_PATH/allsample_merge/POSITION_LIST $DATA_PATH/sampledata/$line/read_info/${line}_data.csv $NCORE
python3 $DATA_PATH/script/add_strand.py $DATA_PATH/sampledata/$line/read_info/${line}_data.csv $DATA_PATH/sampledata/$line/read_info/${line}_*reads_in_cluster.bed > $DATA_PATH/sampledata/$line/read_info/${line}_data_strand.csv
cp $DATA_PATH/sampledata/$line/read_info/${line}_data_strand.csv $DATA_PATH/check_seq/master_file/sampledata/.
echo "=== process6: ${line}の照合完了 ==="
done < $SAMPLE
echo "=== process6: 全サンプルの照合完了 ==="