#!/bin/bash
#bwaインデックスファイル付きのサブジェクトファイル（ERVリードを判別したいデータベースのマルチファスタファイル）を$DATA_PATH/check_seq/bwa/subject直下においておくこと

DATA_PATH=$1
BWA_SUBJECT=$2
NCORE=$3

date
echo "=== process4: fastqのマージとマッピング中 ==="
cat $DATA_PATH/check_seq/fq/sampledata/*.fq > $DATA_PATH/check_seq/fq/allsample_ERVread.fq
bwa mem -t $NCORE -o $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs.bam $BWA_SUBJECT $DATA_PATH/check_seq/fq/allsample_ERVread.fq
samtools view -@ $NCORE -F 2048 -o $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs_F2048.bam $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs.bam
samtools view $DATA_PATH/check_seq/bwa/result/allsample_to_DfamERV_ERVcons_RFs_F2048.bam | cut -f 1-6 - | sort -k1,1 - > $DATA_PATH/check_seq/bwa/result/mapping_info_sort
echo "=== process4: fastqのマージとマッピング完了 ==="
date
echo "=== process5: insertポジションのリストを作成中 ==="
cat $DATA_PATH/allsample_merge/sampledata/*.bed > $DATA_PATH/allsample_merge/allsample_merged.bed
sort -V -k1,1 -k2,2 $DATA_PATH/allsample_merge/allsample_merged.bed | bedtools merge -i - -d 50 > $DATA_PATH/allsample_merge/allsample_merged_d50.bed
awk 'BEGIN {OFS="\t"} {print "POS" NR, $1, $2, $3}' $DATA_PATH/allsample_merge/allsample_merged_d50.bed > $DATA_PATH/allsample_merge/POSITION_LIST
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "POS" NR}' $DATA_PATH/allsample_merge/allsample_merged_d50.bed > $DATA_PATH/allsample_merge/pos.bed
date
echo "=== process5: insertポジションのリスト作成完了 ==="
