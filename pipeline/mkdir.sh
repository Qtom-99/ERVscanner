#!/bin/bash

SAMPLE=$1
DATA_PATH=$2

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
