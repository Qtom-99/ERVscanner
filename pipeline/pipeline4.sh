#!/bin/bash

#あらかじめ以下のようなタブ区切りのヘッダー付き辞書ファイルを$DATA_PATH/dfam_infoの直下にDfam_ERV_infoという名前でおいておくこと（https://dfam.org/browseのフィルタリング後の内容をそのままコピペするとよい）
#Accession	Name	Classification	Clades	Description	Length
#DF000001893	LTRIS2	ERV1	Mus	Mouse subfamily of LTR retrotransposons	564
#第四引数は、ほしいERVclassや名前のリスト（辞書の３列目のグループ名のこと　ERV1、HERVKなど）

DATA_PATH=$1
threshold=$2
ALLSAMPLELIST=$3
ACCESSIONLIST=$4

date
echo "=== process7: 挿入配列の推定によるフィルタリング開始  ==="
cat $DATA_PATH/check_seq/master_file/sampledata/*strand.csv > $DATA_PATH/check_seq/master_file/allsample_data_merge.csv
sort -k5,5 $DATA_PATH/check_seq/master_file/allsample_data_merge.csv > $DATA_PATH/check_seq/master_file/allsample_data_merge_sort.csv
paste $DATA_PATH/check_seq/master_file/allsample_data_merge_sort.csv $DATA_PATH/check_seq/bwa/result/mapping_info_sort > $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master.txt
sort -V -k1,1 -k11,11 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master.txt > $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt
cut -f 1,14 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt > $DATA_PATH/check_seq/master_file/allsample_pos_subject
awk -F'\t' '{split($2, a, "::"); $2 = a[1]; print $1 "\t" $2}' $DATA_PATH/check_seq/master_file/allsample_pos_subject > $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr
python3 $DATA_PATH/script/ERVname_to_class.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr $DATA_PATH/dfam_info/Dfam_ERV_info $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_sed
python3 $DATA_PATH/script/DFaccession_to_class.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_sed $DATA_PATH/dfam_info/Dfam_ERV_info $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_class
python3 $DATA_PATH/script/major_class.py $DATA_PATH/check_seq/master_file/allsample_pos_subject_nochr_class $DATA_PATH/check_seq/master_file/pos_allclass_report $DATA_PATH/check_seq/master_file/pos_majorclass_report
awk -F'\t' -v threshold="$threshold" '$3 >= threshold' "$DATA_PATH/check_seq/master_file/pos_majorclass_report" > "$DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up"
python3 $DATA_PATH/script/extract_wanted_accession.py $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up $ACCESSIONLIST $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted
cut -f 1 $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted | sort -V > $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted_sort
cut -f 1 $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted | sort -V > $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted_sort
sort -V -k1,1 $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted > $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted_sort_ID
cut -f 1,10 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt > $DATA_PATH/check_seq/master_file/pos_sample
sort -V -k1,1 -k2,2 $DATA_PATH/check_seq/master_file/pos_sample | uniq > $DATA_PATH/check_seq/master_file/pos_sample_uniq
cut -f 1,10,11,13 $DATA_PATH/check_seq/master_file/allsample_ERVinfo_master_sort.txt > $DATA_PATH/check_seq/master_file/allsample_POSID_strand_seq
python3 $DATA_PATH/script/get_strand_majoority.py $DATA_PATH/check_seq/master_file/allsample_POSID_strand_seq $DATA_PATH/check_seq/master_file/allsample_POSID_strand_pass $DATA_PATH/check_seq/master_file/allsample_POSID_strand_filtered
python3 $DATA_PATH/script/make01table.py $DATA_PATH/check_seq/master_file/pos_sample_uniq $ALLSAMPLELIST $DATA_PATH/01table/allsample_01table.tsv
python3 $DATA_PATH/script/extract_major_pos.py $DATA_PATH/01table/allsample_01table.tsv $DATA_PATH/check_seq/master_file/pos_majorclass_${threshold//./}up_wanted_sort $DATA_PATH/01table/allsample_01table_${threshold//./}up.tsv
date
echo "=== process7: フィルタリング完了 ${threshold}以上の一致率の挿入のみを抽出しました  ==="