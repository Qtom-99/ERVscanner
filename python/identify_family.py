#!/usr/bin/env python3
import sys
import argparse
import pysam
from collections import defaultdict
import pandas as pd
import numpy as np

# your columns (0-based) based on the example
POS_COL = 0
QNAME_COL = 4
STRAND_COL = 10

def split_auto(line: str):
    # accept TSV or CSV
    if "\t" in line:
        return line.split("\t")
    return line.split(",")

def load_qname_to_pos_lr(path: str):
    q2plr = defaultdict(list)
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = split_auto(line)
            if len(cols) <= max(POS_COL, QNAME_COL, STRAND_COL):
                continue

            pos_id = cols[POS_COL]
            qname = cols[QNAME_COL].split("/")[0]
            strand = cols[STRAND_COL].strip()

            if strand == "+":
                lr = "L"
            elif strand == "-":
                lr = "R"
            else:
                continue

            q2plr[qname].append((pos_id, lr))
    return q2plr

def upd_minmax(mm, a, b):
    # mm is [min, max]; allow init by None
    if mm[0] is None or a < mm[0]:
        mm[0] = a
    if mm[1] is None or b > mm[1]:
        mm[1] = b

def compute_repeat_lr_stats(
    bam_path: str,
    strand_file: str,
    allow_secondary: bool = False,
    pair_mode: str = "any",  # "any" or "read1"
):
    """
    Stream BAM (qname-sorted) and aggregate min/max ranges + counts per (POS, repeat).

    Returns
    -------
    stats : dict
        key: (pos_id, repeat)
        value: dict with counts and min/max ranges
    """
    q2plr = load_qname_to_pos_lr(strand_file)

    stats = defaultdict(lambda: {
        "L_plus_minmax":  [None, None],   # [min_start, max_end)
        "L_minus_minmax": [None, None],
        "R_plus_minmax":  [None, None],
        "R_minus_minmax": [None, None],
        "n_L_plus": 0,
        "n_L_minus": 0,
        "n_R_plus": 0,
        "n_R_minus": 0,
        "ref_len": None,  # repeat length (from BAM header)
    })

    bam = pysam.AlignmentFile(bam_path, "rb", require_index=False)

    for r in bam:  # index-free streaming
        if (not allow_secondary) and (r.is_secondary or r.is_supplementary):
            continue
        if pair_mode == "read1" and (not r.is_read1):
            continue
        if r.is_unmapped:
            continue

        repeat = r.reference_name
        if repeat is None:
            continue

        qname = r.query_name
        if qname not in q2plr:
            continue

        rs = r.reference_start
        re = r.reference_end
        if rs is None or re is None:
            continue

        # repeat length from BAM header
        try:
            ref_len = bam.get_reference_length(repeat)
        except Exception:
            ref_len = None

        # ★重要：pos_id, lr ごとに st を取って、その場でカウント/更新する
        for pos_id, lr in q2plr[qname]:
            key = (pos_id, repeat)
            st = stats[key]

            if st["ref_len"] is None:
                st["ref_len"] = ref_len

            if lr == "L":
                if r.is_reverse:
                    upd_minmax(st["L_minus_minmax"], rs, re)
                    st["n_L_minus"] += 1
                else:
                    upd_minmax(st["L_plus_minmax"], rs, re)
                    st["n_L_plus"] += 1

            elif lr == "R":
                if r.is_reverse:
                    upd_minmax(st["R_minus_minmax"], rs, re)
                    st["n_R_minus"] += 1
                else:
                    upd_minmax(st["R_plus_minmax"], rs, re)
                    st["n_R_plus"] += 1

            else:
                # lr が想定外なら無視（必要なら raise にしてもOK）
                continue

    bam.close()
    return stats


def write_repeat_lr_stats_tsv(stats, out_tsv: str):
    header = (
        "POS\trepeat\trep_len\t"
        "n_L_plus\tn_L_minus\tn_R_plus\tn_R_minus\t"
        "L_plus_min\tL_plus_max\tL_minus_min\tL_minus_max\t"
        "R_plus_min\tR_plus_max\tR_minus_min\tR_minus_max\n"
    )

    with open(out_tsv, "w") as out:
        out.write(header)

        for (pos_id, repeat), st in sorted(stats.items()):
            rep_len = st["ref_len"] if st["ref_len"] is not None else "NA"

            L_plus_min,  L_plus_max  = st["L_plus_minmax"]
            L_minus_min, L_minus_max = st["L_minus_minmax"]
            R_plus_min,  R_plus_max  = st["R_plus_minmax"]
            R_minus_min, R_minus_max = st["R_minus_minmax"]

            out.write(
                f"{pos_id}\t{repeat}\t{rep_len}\t"
                f"{st['n_L_plus']}\t{st['n_L_minus']}\t{st['n_R_plus']}\t{st['n_R_minus']}\t"
                f"{L_plus_min}\t{L_plus_max}\t{L_minus_min}\t{L_minus_max}\t"
                f"{R_plus_min}\t{R_plus_max}\t{R_minus_min}\t{R_minus_max}\n"
            )


def run_pipeline(bam_path: str, strand_file: str, out_tsv: str,
                 allow_secondary: bool = False, pair_mode: str = "any"):
    stats = compute_repeat_lr_stats(
        bam_path=bam_path,
        strand_file=strand_file,
        allow_secondary=allow_secondary,
        pair_mode=pair_mode,
    )
    write_repeat_lr_stats_tsv(stats, out_tsv)
    return out_tsv

def read_repeat_lr_stats_tsv(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", dtype={"POS": str, "repeat": str})
    # 数値列を安全に数値化（None/NA などは NaN）
    num_cols = [
        "rep_len",
        "n_L_plus","n_L_minus","n_R_plus","n_R_minus",
        "L_plus_min","L_plus_max","L_minus_min","L_minus_max",
        "R_plus_min","R_plus_max","R_minus_min","R_minus_max",
    ]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # count は NaN があり得ない想定だが、一応 0 埋め
    for c in ["n_L_plus","n_L_minus","n_R_plus","n_R_minus"]:
        if c in df.columns:
            df[c] = df[c].fillna(0).astype(int)
    return df

def call_insertion_repeat(
    df: pd.DataFrame,
    min_support: int = 4,      # <=3 は挿入なし
    max_span: int = 1000,      # 採用ペアのL or Rのレンジ幅が >=1000 なら挿入なし
) -> pd.DataFrame:
    """
    For each POS, select repeat with maximum support:
      s1 = min(n_L_plus,  n_R_minus)
      s2 = min(n_L_minus, n_R_plus)
      support = max(s1, s2)

    Adopted pair:
      - if s1 > s2 => L+_R-
      - if s2 > s1 => L-_R+
      - if tie => choose by (n_L_plus+n_R_minus) vs (n_L_minus+n_R_plus), tie -> L+_R-

    Then:
      - strand: L+_R- => '-', L-_R+ => '+'
      - spans: adopt_L_span, adopt_R_span; if either >= max_span => NA
      - ins_len computed only from adopted pair ranges:
            all_min = min(adopt_L_min, adopt_R_min)
            all_max = max(adopt_L_max, adopt_R_max)
            ins_len = all_max - all_min  (only if both sides exist)
    """
    df = df.copy()

    # support
    df["support_Lplus_Rminus"] = np.minimum(df["n_L_plus"], df["n_R_minus"])
    df["support_Lminus_Rplus"] = np.minimum(df["n_L_minus"], df["n_R_plus"])
    df["support"] = np.maximum(df["support_Lplus_Rminus"], df["support_Lminus_Rplus"])

    # tie-break helper for adopted pair
    df["pair_score_1"] = df["n_L_plus"] + df["n_R_minus"]   # for L+_R-
    df["pair_score_2"] = df["n_L_minus"] + df["n_R_plus"]   # for L-_R+

    # tie-break among repeats
    df["total_reads"] = df["n_L_plus"] + df["n_L_minus"] + df["n_R_plus"] + df["n_R_minus"]
    df["_rep_len_sort"] = df["rep_len"].fillna(-1)

    # POSごとにベストrepeat
    best = (
        df.sort_values(
            by=["POS", "support", "total_reads", "_rep_len_sort", "repeat"],
            ascending=[True, False, False, False, True],
            kind="mergesort",
        )
        .groupby("POS", as_index=False)
        .head(1)
        .copy()
    )

    # adopted pair（支持度ロジックで決定）
    best["adopted_pair"] = np.where(
        best["support_Lplus_Rminus"] > best["support_Lminus_Rplus"], "L+_R-",
        np.where(
            best["support_Lplus_Rminus"] < best["support_Lminus_Rplus"], "L-_R+",
            np.where(best["pair_score_1"] >= best["pair_score_2"], "L+_R-", "L-_R+")
        )
    )

    # adopted L/R min/max（採用ペアに対応）
    best["adopt_L_min"] = np.where(best["adopted_pair"] == "L+_R-", best["L_plus_min"],  best["L_minus_min"])
    best["adopt_L_max"] = np.where(best["adopted_pair"] == "L+_R-", best["L_plus_max"],  best["L_minus_max"])
    best["adopt_R_min"] = np.where(best["adopted_pair"] == "L+_R-", best["R_minus_min"], best["R_plus_min"])
    best["adopt_R_max"] = np.where(best["adopted_pair"] == "L+_R-", best["R_minus_max"], best["R_plus_max"])

    best["adopt_L_span"] = best["adopt_L_max"] - best["adopt_L_min"]
    best["adopt_R_span"] = best["adopt_R_max"] - best["adopt_R_min"]

    # 挿入向き（採用ペア由来）
    best["ins_strand"] = np.where(best["adopted_pair"] == "L+_R-", "-", "+")

    # 採用ペアのみで insertion length
    best["pair_all_min"] = best[["adopt_L_min", "adopt_R_min"]].min(axis=1, skipna=True)
    best["pair_all_max"] = best[["adopt_L_max", "adopt_R_max"]].max(axis=1, skipna=True)
    best["ins_len"] = np.where(
        best["pair_all_min"].isna() | best["pair_all_max"].isna(),
        np.nan,
        best["pair_all_max"] - best["pair_all_min"]
    )

    # 判定フラグ
    best["pass_support"] = best["support"] >= min_support
    best["pass_span"] = ~(
        (best["adopt_L_span"] >= max_span) |
        (best["adopt_R_span"] >= max_span)
    )
    best["pass_len"] = best["ins_len"] >= 200

    # call（挿入なしなら NA / NaN）
    ok = best["pass_support"] & best["pass_span"] & best["pass_len"]

    best["called_repeat"] = np.where(ok, best["repeat"], "NA")
    best["called_support"] = best["support"]
    best["called_pair"] = np.where(ok, best["adopted_pair"], "NA")
    best["called_strand"] = np.where(ok, best["ins_strand"], "NA")

    # 付随数値（NGなら NaN）
    for c in ["adopt_L_min","adopt_L_max","adopt_R_min","adopt_R_max",
              "adopt_L_span","adopt_R_span","pair_all_min","pair_all_max","ins_len"]:
        best[c] = np.where(ok, best[c], np.nan)

    out = best[[
        "POS",
        "called_repeat", "called_support", "called_pair", "called_strand",
        "adopt_L_min","adopt_L_max","adopt_R_min","adopt_R_max",
        "adopt_L_span","adopt_R_span",
        "pair_all_min","pair_all_max","ins_len",
        # 監査用
        "repeat","rep_len",
        "n_L_plus","n_L_minus","n_R_plus","n_R_minus",
        "support_Lplus_Rminus","support_Lminus_Rplus","support","total_reads",
    ]].rename(columns={"repeat": "best_repeat", "rep_len": "best_rep_len"})

    return out

def run_calling_from_tsv(in_tsv: str, out_call_tsv: str, min_support: int = 4, max_span: int = 1000):
    df = read_repeat_lr_stats_tsv(in_tsv)
    called = call_insertion_repeat(df, min_support=min_support, max_span=max_span)
    called = called[called["called_repeat"] != 'NA']
    
    called.to_csv(out_call_tsv, sep="\t", index=False)
    called_pos = called["POS"]
    base = out_call_tsv.removesuffix(".tsv")
    called_pos.to_csv(base+"_pos.txt", sep="\t", index=False, header=False)
    
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam", help="qname-sorted BAM mapped to repeats (reference_name=repeat)")
    ap.add_argument("strand_file", help="strand.tsv or .csv (POS, qname, strand)")
    ap.add_argument("out_tsv", help="output TSV")
    ap.add_argument("--allow-secondary", action="store_true", default=True,
                    help="include secondary/supplementary (default: include)")
    ap.add_argument("--pair-mode", choices=["any", "read1"], default="any",
                    help="any: count any alignment (default). read1: only read1 (if flags are reliable).")
    ap.add_argument("--min-support", type=int, default=4,
                help="minimum support to call insertion (<=3 => NA). default: 4")
    ap.add_argument("--max-span", type=int, default=1000,
                help="if adopted L or R span >= this, call NA. default: 1000")
    ap.add_argument("--out-call", default=None,
                help="output TSV for called insertion per POS (default: out_tsv + '.call.tsv')")
    args = ap.parse_args()

    raw_tsv = run_pipeline(
        bam_path=args.bam,
        strand_file=args.strand_file,
        out_tsv=args.out_tsv,
        allow_secondary=args.allow_secondary,
        pair_mode=args.pair_mode,
    )
    base = raw_tsv.removesuffix(".tsv")
    out_call = args.out_call if args.out_call is not None else (base + ".call.tsv")
    run_calling_from_tsv(raw_tsv, out_call_tsv=out_call, min_support=args.min_support, max_span=args.max_span)


if __name__ == "__main__":
    main()