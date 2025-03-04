#!/usr/bin/env python3
import argparse
import math
import os
import pandas as pd

def read_breakpoints(breakpoint_file, bed_file):
    """
    ファイル１（breakpoints.tsv）とファイル２（regions.bed）を読み込み、
    POSIDごとに (chrom, pos, qual) の情報を辞書として返す。
      - ファイル１にあればその情報 (QUAL="high")
      - ない場合は、BEDの start/end の平均（四捨五入）を採用し QUAL="low"
    """
    pos_breakpoints = {}
    # ファイル１の読み込み
    with open(breakpoint_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, pos_str, posid = parts[0], parts[1], parts[2]
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            pos_breakpoints[posid] = (chrom, pos, "high")
    
    # ファイル２（BED）の読み込み
    with open(bed_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start_str, end_str, posid = parts[0], parts[1], parts[2], parts[3]
            if posid not in pos_breakpoints:
                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    continue
                bp = int(round((start + end) / 2))
                pos_breakpoints[posid] = (chrom, bp, "low")
    return pos_breakpoints

def read_mapping(file_path, key_index=0, value_index=1):
    """
    汎用関数：タブ区切りファイルから key_index 列をキー、value_index 列を値とした dict を返す。
    """
    mapping = {}
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= value_index:
                continue
            key = parts[key_index]
            value = parts[value_index]
            mapping[key] = value
    return mapping

def load_sample_coverage(sample_ids, cov_dir):
    """
    ファイル８のディレクトリから、各サンプルのカバードリード情報を読み込み、
    sample_cov[sample_id] = { POSID: CR } の辞書を返す。
    各ファイルは "<cov_dir>/<SAMPLEID>_n_of_covreads_noclip.txt" という名前とする。
    ファイルが存在しない場合は空の辞書を設定。
    """
    sample_cov = {}
    for sample in sample_ids:
        filename = os.path.join(cov_dir, f"{sample}_n_of_covreads_noclip.txt")
        mapping = {}
        if os.path.exists(filename):
            with open(filename, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 5:
                        continue
                    posid_line = parts[3]
                    cr_val = parts[4]
                    mapping[posid_line] = cr_val
        # 存在しない場合は空の辞書
        sample_cov[sample] = mapping
    return sample_cov

def process_single_posid(posid, genotype_values, sample_ids, sample_cov, pos_breakpoints, pos_ervclass, pos_ervfamily, pos_direction):
    """
    各POSIDの処理を行い、VCFレコード1行の文字列を返す。
    genotype_values: ファイル３の該当POSID列の各サンプルの値（文字列のリスト）
    sample_ids: サンプルIDのリスト（ファイル３またはファイル４の順序）
    sample_cov: 各サンプルのカバードリード情報の辞書（sample_id -> {POSID: CR}）
    """
    ac = 0  # ALTアレル出現数
    an = 0  # コールされた総アレル数（ジェノタイプ不能は除く）
    sample_fields = []  # 各サンプルの "GT:CR" フィールド
    for i, val in enumerate(genotype_values):
        try:
            v = int(val)
        except ValueError:
            v = 0
        if v == 0:
            gt = "0/0"
            an += 2
            cr = "."
        elif v == 1:
            gt = "0/1"
            ac += 1
            an += 2
            # サンプルのcoverage情報から、このPOSIDのCRを取得
            cr = sample_cov.get(sample_ids[i], {}).get(posid, ".")
        elif v == 2:
            gt = "1/1"
            ac += 2
            an += 2
            cr = sample_cov.get(sample_ids[i], {}).get(posid, ".")
        elif v == 3:
            gt = "1/."
            # 挿入はあるがジェノタイピング不能の場合、CRは情報なし
            cr = "."
        else:
            gt = "./."
            cr = "."
        sample_fields.append(f"{gt}:{cr}")
    af = (ac / an) if an > 0 else 0

    if posid not in pos_breakpoints:
        print(f"Warning: POSID {posid} がブレークポイント情報に存在しません。スキップします。")
        return None
    chrom, pos, qual = pos_breakpoints[posid]
    ec = pos_ervclass.get(posid, ".")
    ef = pos_ervfamily.get(posid, ".")
    dr = pos_direction.get(posid, ".")
    info = f"AC={ac};AF={af:.3f};AN={an};DP=0;EC={ec};EF={ef};DR={dr}"
    # FORMAT は GT:CR
    fields = [chrom, str(pos), posid, ".", "INS_seq", qual, ".", info, "GT:CR"]
    fields.extend(sample_fields)
    return "\t".join(fields)

def main():
    parser = argparse.ArgumentParser(
        description="各種入力ファイルから情報を統合し、シングルコアでVCFファイルを作成します。\n"
                    "ファイル３の横軸からレコードすべきPOSID一覧を取得し、\n"
                    "ファイル１にあればその情報を、なければファイル２からブレークポイントを求めます。\n"
                    "さらに、ファイル８の各サンプルのカバードリード情報をGTの横にCRとして追加します。"
    )
    parser.add_argument("--file1", required=True, help="ファイル１: ブレークポイント情報 (TSV)")
    parser.add_argument("--file2", required=True, help="ファイル２: 領域情報BED (TSV)")
    parser.add_argument("--file3", required=True, help="ファイル３: サンプル×POSID のジェノタイプTSV")
    parser.add_argument("--file4", help="ファイル４: サンプルIDリスト（オプション）")
    parser.add_argument("--file5", required=True, help="ファイル５: POSID と ERVCLASS (TSV)")
    parser.add_argument("--file6", required=True, help="ファイル６: POSID と ERVFAMILY (TSV)")
    parser.add_argument("--file7", required=True, help="ファイル７: POSID と 方向 (TSV)")
    parser.add_argument("--file8", required=True, help="ファイル８: サンプルごとのカバードリード情報が格納されたディレクトリ")
    parser.add_argument("--output", required=True, help="出力VCFファイル")
    args = parser.parse_args()

    # --- ブレークポイント情報の取得 ---
    pos_breakpoints = read_breakpoints(args.file1, args.file2)

    # --- マッピング情報の読み込み ---
    pos_ervclass  = read_mapping(args.file5, key_index=0, value_index=1)  # EC
    pos_ervfamily = read_mapping(args.file6, key_index=0, value_index=1)  # EF
    pos_direction = read_mapping(args.file7, key_index=0, value_index=1)  # DR

    # --- ファイル３のジェノタイプ行列の読み込み ---
    df_geno = pd.read_csv(args.file3, sep="\t", dtype=str)
    sample_ids_from_file3 = df_geno.iloc[:, 0].tolist()

    # --- ファイル４が指定されていれば、サンプル順をそれに合わせる ---
    if args.file4:
        with open(args.file4, "r", encoding="utf-8") as f:
            sample_order = [line.strip() for line in f if line.strip()]
        df_geno = df_geno.set_index(df_geno.columns[0])
        df_geno = df_geno.reindex(sample_order).dropna().reset_index()
        sample_ids = sample_order
    else:
        sample_ids = sample_ids_from_file3

    # --- ファイル８から各サンプルのカバードリード情報を読み込む ---
    sample_cov = load_sample_coverage(sample_ids, args.file8)

    # --- ファイル３の横軸からレコードすべきPOSID一覧を取得 ---
    posid_list = list(df_geno.columns[1:])

    # --- 各POSIDごとに、各サンプルのジェノタイプ値を辞書化 ---
    genotype_by_pos = {}
    for posid in posid_list:
        genotype_by_pos[posid] = df_geno[posid].tolist()

    # --- 各POSIDごとにVCFレコードを生成（シングルコア処理） ---
    results = []
    for posid in posid_list:
        record = process_single_posid(posid, genotype_by_pos[posid], sample_ids, sample_cov,
                                      pos_breakpoints, pos_ervclass, pos_ervfamily, pos_direction)
        if record is not None:
            results.append(record)

    # --- VCFファイルの出力 ---
    with open(args.output, "w", encoding="utf-8") as outvcf:
        # VCFヘッダーの出力
        outvcf.write("##fileformat=VCFv4.2\n")
        outvcf.write("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"ALT allele count\">\n")
        outvcf.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"ALT allele frequency\">\n")
        outvcf.write("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n")
        outvcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n")
        outvcf.write("##INFO=<ID=EC,Number=1,Type=String,Description=\"ERVCLASS\">\n")
        outvcf.write("##INFO=<ID=EF,Number=1,Type=String,Description=\"ERVFAMILY\">\n")
        outvcf.write("##INFO=<ID=DR,Number=1,Type=String,Description=\"Insertion direction\">\n")
        outvcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        outvcf.write("##FORMAT=<ID=CR,Number=1,Type=Integer,Description=\"Covered read count\">\n")
        header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids) + "\n"
        outvcf.write(header_line)
        for rec in results:
            outvcf.write(rec + "\n")

    print("VCFの作成が完了しました。出力ファイル:", args.output)

if __name__ == "__main__":
    main()
