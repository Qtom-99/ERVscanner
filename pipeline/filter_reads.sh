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

REF_GENOME="$DATA_PATH/reference/reference.fasta"
QUERY_BED="dfam_info/target.bed"
ALT_CHR_LIST="$DATA_PATH/reference/alt_chr_list"


# reading confing file
eval $(awk -F: '
    # コメント・空行をスキップ
    /^[[:space:]]*$/ {next}      # skip empty line
    /^[[:space:]]*#/ {next}      # ignore comment line

    # remove spaces
    {gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)}

    # exporting parameters
    {print $1"=\"" $2 "\""}
' "$CONFIG_FILE")

while read line
do
date
echo "=== process1: SAMPLE:${line} estracting ERV regions ==="
samtools view -@ "$NCORE" -T "$REF_GENOME" -L $DATA_PATH/$QUERY_BED -P -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap.bam -O BAM "$INPUT_PATH"/${line}."$INPUT_TYPE"
echo "=== process1: SAMPLE:${line} estracting ERV regions: finished ==="
date
echo "=== process2: SAMPLE:${line} insertion estimation ==="
samtools view -@ $NCORE -b -F 256 $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap.bam |samtools view -@ $NCORE -F 2048 -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048.bam -
samtools sort -@ $NCORE -n -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048.bam
bedtools pairtobed -type xor -abam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn.bam -b $DATA_PATH/$QUERY_BED 2>/dev/null > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam
bedtools bamtobed -i $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bed
awk -F'\t' '$4 ~ "/1$"' $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bed > $DATA_PATH/sampledata/${line}/${line}_read1.bed
awk -F'\t' '$4 ~ "/2$"' $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bed > $DATA_PATH/sampledata/${line}/${line}_read2.bed
python3 $PY_PATH/bed_curate.py $DATA_PATH/sampledata/${line}/${line}_read1.bed $DATA_PATH/sampledata/${line}/${line}_read2.bed $DATA_PATH/sampledata/${line}/${line}_pair.bed
python3 $PY_PATH/other_chr.py $DATA_PATH/sampledata/${line}/${line}_pair.bed $DATA_PATH/sampledata/${line}/${line}_other_chr.bed
python3 $PY_PATH/samechr_split.py $DATA_PATH/sampledata/${line}/${line}_pair.bed $DATA_PATH/sampledata/${line}/${line}_samechr_split.bed
python3 $PY_PATH/extract_inproper_pair.py -i $DATA_PATH/sampledata/${line}/${line}_pair.bed -o $DATA_PATH/sampledata/${line}/${line}_inproper.bed
cat $DATA_PATH/sampledata/${line}/${line}_other_chr.bed $DATA_PATH/sampledata/${line}/${line}_inproper.bed $DATA_PATH/sampledata/${line}/${line}_samechr_split.bed > $DATA_PATH/sampledata/${line}/${line}_mateunmapped.bed
cut -f 4 $DATA_PATH/sampledata/${line}/${line}_mateunmapped.bed > $DATA_PATH/sampledata/${line}/${line}_readname.txt
python3 $PY_PATH/remove_endname.py $DATA_PATH/sampledata/${line}/${line}_readname.txt $DATA_PATH/sampledata/${line}/${line}_readname_nonread.txt
sort $DATA_PATH/sampledata/${line}/${line}_readname_nonread.txt | uniq - > $DATA_PATH/sampledata/${line}/${line}_readname_nonread_uniq.txt
samtools view -@ $NCORE -N $DATA_PATH/sampledata/${line}/${line}_readname_nonread_uniq.txt -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam
bedtools bamtobed -i $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bam > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bed
bedtools intersect -v -abam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bam -b $DATA_PATH/$QUERY_BED > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8.bam
samtools sort -@ $NCORE -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8.bam
samtools view -q ${MAPQ_THRESHOLD} -o $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort.bam
python3 $PY_PATH/extract_altxahit.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}.bam $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}_noalt.bam $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}_remove.bam $ALT_CHR_LIST
bedtools bamtobed -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}_noalt.bam > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}_noalt.bed
bedtools sort -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}_noalt.bed | bedtools cluster -s -d 200 -i - > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200.bed
python3 $PY_PATH/bed_cluster.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads.bed $CLUSTER_THRESHOLD
bedtools sort -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads.bed | bedtools merge -i - -d 200 -s -c 6 -o distinct > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged.bed
sort -V -k1,1 -k2,2 $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged.bed > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort.bed
bedtools cluster -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort.bed -d 200 > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered.bed
python3 $PY_PATH/bed_cluster_plusminus.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus.bed
python3 $PY_PATH/bed_cluster_plusminus_jun.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun.bed
python3 $PY_PATH/bed_merge2.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun_merged2.bed
cp $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun_merged2.bed $DATA_PATH/allsample_merge/sampledata/.
echo "=== process2: SAMPLE:${line} insertion estimation: finished ==="
date
echo "=== process3: SAMPLE:${line} generating fastq files of ERV regions ==="
bedtools intersect -a $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads.bed -b $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${MAPQ_THRESHOLD}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun.bed > $DATA_PATH/sampledata/${line}/read_info/${line}_${CLUSTER_THRESHOLD}reads_in_cluster.bed
cut -f 4 $DATA_PATH/sampledata/${line}/read_info/${line}_${CLUSTER_THRESHOLD}reads_in_cluster.bed > $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname.txt
python3 $PY_PATH/remove_endname.py $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname.txt $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend.txt
uniq $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend.txt > $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend_uniq
bedtools intersect -a $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam -b $DATA_PATH/$QUERY_BED -v > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition.bam
samtools view -@ $NCORE -F 256 -b $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition.bam | samtools view -@ $NCORE -b -F 2048 -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_F256_F2048.bam -
samtools view -@ $NCORE -N $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend_uniq -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.bam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_F256_F2048.bam
samtools view -@ $NCORE -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.sam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.bam
cut -f 1,3,4,7,8 $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.sam > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut.sam
awk 'BEGIN{FS=OFS="\t"} $4=="="{ $4=$2 }1' $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut.sam > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut_awk.sam
sort -V -k2,2 -k3,3 $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut_awk.sam > $DATA_PATH/sampledata/${line}/read_info/${line}
bedtools intersect -a $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam -b $DATA_PATH/$QUERY_BED -wa > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread.bam
samtools view -@ $NCORE -F 256 -b $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread.bam | samtools view -@ $NCORE -b -F 2048 -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_F256_F2048.bam -
samtools view -@ $NCORE -N $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend_uniq -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true.bam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_F256_F2048.bam
samtools sort -@ $NCORE -n -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.bam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true.bam
bedtools bamtofastq -i $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.bam -fq $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.fq
cp $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.fq $DATA_PATH/check_seq/fq/sampledata/.
rm $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048.bam
rm $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn.bam
rm $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bed
rm $DATA_PATH/sampledata/${line}/${line}_pair.bed
rm $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_F256_F2048.bam
rm $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_F256_F2048.bam
echo "=== process3: SAMPLE:${line} generating fastq files of ERV regions: finished ==="
done < $SAMPLE
