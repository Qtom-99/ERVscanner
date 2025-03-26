#!/bin/bash
# samtools 1.20 or higher is recommended
# for input files, please see README.md file

# default values
QUALITY=30
CLUSTER_THRESHOLD=5
INPUT_TYPE="BAM"

# getting option values
while getopts "i:o:" opt; do
  case $opt in
    s) SAMPLE="$OPTARG" ;;
    i) INPUT_PATH="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    t) INPUT_TYPE="$OPTARG" ;;
    d) DATA_PATH="$OPTARG" ;;
    n) NCORE="$OPTARG" ;;
    b) QUERY_BED_FULLPATH="$OPTARG" ;;
    q) QUALITY="$OPTARG" ;;
    c) CLUSTER_THRESHOLD="$OPTARG" ;;
    a) ALT_CHR_LIST="$OPTARG" ;;
    \?) echo "Usage: $0 [-i input] [-o output]" >&2; exit 1 ;;
  esac
done

if [[ -z "$input" ]]; then
  echo "Error: -i is required" >&2
  exit 1
fi

while read line
do
date
echo "=== process1: SAMPLE:${line} ERV領域抽出開始 ==="
samtools view -@ "$NCORE" -T "$REF_GENOME" -L $QUERY_BED_FULLPATH -P -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap.bam -O BAM "$INPUT_PATH"/${line}."$INPUT_TYPE"
echo "=== process1: SAMPLE:${line} ERV領域抽出完了 ==="
date
echo "=== process2: SAMPLE:${line} insertion推定開始 ==="
samtools view -@ $NCORE -b -F 256 $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap.bam |samtools view -@ $NCORE -F 2048 -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048.bam -
samtools sort -@ $NCORE -n -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048.bam
bedtools pairtobed -type xor -abam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn.bam -b $QUERY_BED_FULLPATH 2>/dev/null > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam
bedtools bamtobed -i $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bed
awk -F'\t' '$4 ~ "/1$"' $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bed > $DATA_PATH/sampledata/${line}/${line}_read1.bed
awk -F'\t' '$4 ~ "/2$"' $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bed > $DATA_PATH/sampledata/${line}/${line}_read2.bed
python3 $DATA_PATH/script/bed_curate.py $DATA_PATH/sampledata/${line}/${line}_read1.bed $DATA_PATH/sampledata/${line}/${line}_read2.bed $DATA_PATH/sampledata/${line}/${line}_pair.bed
python3 $DATA_PATH/script/other_chr.py $DATA_PATH/sampledata/${line}/${line}_pair.bed $DATA_PATH/sampledata/${line}/${line}_other_chr.bed
python3 $DATA_PATH/script/samechr_split.py $DATA_PATH/sampledata/${line}/${line}_pair.bed $DATA_PATH/sampledata/${line}/${line}_samechr_split.bed
python3 $DATA_PATH/script/extract_inproper_pair.py -i $DATA_PATH/sampledata/${line}/${line}_pair.bed -o $DATA_PATH/sampledata/${line}/${line}_inproper.bed
cat $DATA_PATH/sampledata/${line}/${line}_other_chr.bed $DATA_PATH/sampledata/${line}/${line}_inproper.bed $DATA_PATH/sampledata/${line}/${line}_samechr_split.bed > $DATA_PATH/sampledata/${line}/${line}_mateunmapped.bed
cut -f 4 $DATA_PATH/sampledata/${line}/${line}_mateunmapped.bed > $DATA_PATH/sampledata/${line}/${line}_readname.txt
python3 $DATA_PATH/script/remove_endname.py $DATA_PATH/sampledata/${line}/${line}_readname.txt $DATA_PATH/sampledata/${line}/${line}_readname_nonread.txt
sort $DATA_PATH/sampledata/${line}/${line}_readname_nonread.txt | uniq - > $DATA_PATH/sampledata/${line}/${line}_readname_nonread_uniq.txt
samtools view -@ $NCORE -N $DATA_PATH/sampledata/${line}/${line}_readname_nonread_uniq.txt -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam
bedtools bamtobed -i $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bam > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bed
bedtools intersect -v -abam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor_discordant.bam -b $QUERY_BED_FULLPATH > $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8.bam
samtools sort -@ $NCORE -o $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8.bam
samtools view -q ${QUALITY} -o $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}.bam $DATA_PATH/sampledata/${line}/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort.bam
python3 $DATA_PATH/script/extract_altxahit.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}.bam $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${QUALITY}_noalt.bam $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${QUALITY}_remove.bam $ALT_CHR_LIST
bedtools bamtobed -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${QUALITY}_noalt.bam > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${QUALITY}_noalt.bed
bedtools sort -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${QUALITY}_noalt.bed | bedtools cluster -s -d 200 -i - > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${QUALITY}clustered200.bed
python3 $DATA_PATH/script/bed_cluster.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_both_discordant_onlyf8_sort_q${QUALITY}clustered200.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads.bed $CLUSTER_THRESHOLD
bedtools sort -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads.bed | bedtools merge -i - -d 200 -s -c 6 -o distinct > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged.bed
sort -V -k1,1 -k2,2 $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged.bed > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort.bed
bedtools cluster -i $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort.bed -d 200 > $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered.bed
python3 $DATA_PATH/script/bed_cluster_plusminus.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus.bed
python3 $DATA_PATH/script/bed_cluster_plusminus_jun.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun.bed
python3 $DATA_PATH/script/bed_merge2.py $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun.bed $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun_merged2.bed
cp $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun_merged2.bed $DATA_PATH/allsample_merge/sampledata/.
echo "=== process2: SAMPLE:${line} insertion推定完了 ==="
date
echo "=== process3: SAMPLE:${line} ERVreadのfastqファイルを準備中 ==="
bedtools intersect -a $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads.bed -b $DATA_PATH/sampledata/${line}/curated_bam/${line}_dfamallhit_ERV_bed_position_F256_F2048_sortn_xor_discordant_onlyf8_sort_q${QUALITY}clustered200_${CLUSTER_THRESHOLD}reads_merged_sort_clustered_plusminus_jun.bed > $DATA_PATH/sampledata/${line}/read_info/${line}_${CLUSTER_THRESHOLD}reads_in_cluster.bed
cut -f 4 $DATA_PATH/sampledata/${line}/read_info/${line}_${CLUSTER_THRESHOLD}reads_in_cluster.bed > $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname.txt
python3 $DATA_PATH/script/remove_endname.py $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname.txt $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend.txt
uniq $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend.txt > $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend_uniq
bedtools intersect -a $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam -b $QUERY_BED_FULLPATH -v > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition.bam
samtools view -@ $NCORE -F 256 -b $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition.bam | samtools view -@ $NCORE -b -F 2048 -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_F256_F2048.bam -
samtools view -@ $NCORE -N $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend_uniq -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.bam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_F256_F2048.bam
samtools view -@ $NCORE -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.sam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.bam
cut -f 1,3,4,7,8 $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true.sam > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut.sam
awk 'BEGIN{FS=OFS="\t"} $4=="="{ $4=$2 }1' $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut.sam > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut_awk.sam
sort -V -k2,2 -k3,3 $DATA_PATH/sampledata/${line}/read_info/${line}_ERVposition_read_true_cut_awk.sam > $DATA_PATH/sampledata/${line}/read_info/${line}
bedtools intersect -a $DATA_PATH/sampledata/${line}/${line}_dfamallhit_overlap_F256_F2048_sortn_xor.bam -b $QUERY_BED_FULLPATH -wa > $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread.bam
samtools view -@ $NCORE -F 256 -b $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread.bam | samtools view -@ $NCORE -b -F 2048 -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_F256_F2048.bam -
samtools view -@ $NCORE -N $DATA_PATH/sampledata/${line}/read_info/${line}_posreadname_nonend_uniq -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true.bam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_F256_F2048.bam
samtools sort -@ $NCORE -n -o $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.bam $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true.bam
bedtools bamtofastq -i $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.bam -fq $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.fq
cp $DATA_PATH/sampledata/${line}/read_info/${line}_ERVread_true_sort.fq $DATA_PATH/check_seq/fq/sampledata/.
echo "=== process3: SAMPLE:${line} RERVreadのfastqファイル作成完了 ==="
done < $SAMPLE
