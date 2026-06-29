########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

#paths
INFILE = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/ABCA4_spliceai_match/abca4_spliceai_variant_event_matches.tsv"
OUTDIR = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/ABCA4_spliceai_match/strong_support_summary"
os.makedirs(OUTDIR, exist_ok=True)

OUT_ALL = os.path.join(OUTDIR, "abca4_strong_support_all.tsv")
OUT_VERY_STRONG = os.path.join(OUTDIR, "abca4_very_strong_support.tsv")
OUT_BEST_PER_PBID = os.path.join(OUTDIR, "abca4_strong_support_best_per_pbid.tsv")
OUT_SUMMARY_BY_VARIANT = os.path.join(OUTDIR, "abca4_strong_support_summary_by_variant.tsv")
OUT_REPORT = os.path.join(OUTDIR, "abca4_strong_support_report.txt")


#define functions
def safe_numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


#Read input
df = pd.read_csv(INFILE, sep="\t", dtype=str)

num_cols = [
    "pb_exon_index", "ref_exon_index",
    "observed_start", "observed_end",
    "observed_donor_pos", "observed_acceptor_pos",
    "reference_donor_pos", "reference_acceptor_pos",
    "variant_POS",
    "DS_AG", "DS_AL", "DS_DG", "DS_DL",
    "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    "pred_AG_pos", "pred_AL_pos", "pred_DG_pos", "pred_DL_pos",
    "maxDS",
    "diff_AG_obs_acceptor", "diff_DG_obs_donor",
    "diff_AL_ref_acceptor", "diff_DL_ref_donor",
    "score_AG_obs_acceptor", "score_DG_obs_donor",
    "score_AL_ref_acceptor", "score_DL_ref_donor",
    "match_score"
]
df = safe_numeric(df, num_cols)

rename_map = {}
for c in df.columns:
    if c.endswith("pred_DL_posmaxDS"):
        pass

#exclude terminal exon categories
exclude_event_types = {"alternative_first_exon", "alternative_last_exon"}
df_main = df.loc[~df["event_type"].isin(exclude_event_types)].copy()

#Define strong support rules
#alt_acceptor
aa_loss = (
    (df_main["event_type"] == "alt_acceptor") &
    (df_main["diff_AL_ref_acceptor"] == 0) &
    (df_main["DS_AL"] >= 0.8)
)
aa_gain = (
    (df_main["event_type"] == "alt_acceptor") &
    (df_main["diff_AG_obs_acceptor"] == 0) &
    (df_main["DS_AG"] >= 0.5)
)
aa_very = aa_loss & aa_gain

#alt_donor
ad_loss = (
    (df_main["event_type"] == "alt_donor") &
    (df_main["diff_DL_ref_donor"] == 0) &
    (df_main["DS_DL"] >= 0.8)
)
ad_gain = (
    (df_main["event_type"] == "alt_donor") &
    (df_main["diff_DG_obs_donor"] == 0) &
    (df_main["DS_DG"] >= 0.5)
)
ad_very = ad_loss & ad_gain

#novel_exon: require both ends supported
ne_acceptor = (
    (df_main["event_type"] == "novel_exon") &
    (df_main["diff_AG_obs_acceptor"] == 0) &
    (df_main["DS_AG"] >= 0.5)
)
ne_donor = (
    (df_main["event_type"] == "novel_exon") &
    (df_main["diff_DG_obs_donor"] == 0) &
    (df_main["DS_DG"] >= 0.5)
)
ne_very = ne_acceptor & ne_donor
ne_strong = ne_very.copy()

#exon_skipping / intron retention
strong_mask = aa_loss | aa_gain | ad_loss | ad_gain | ne_strong
very_strong_mask = aa_very | ad_very | ne_very

strong = df_main.loc[strong_mask].copy()
very_strong = df_main.loc[very_strong_mask].copy()

#Add interpretation columns
def classify_row(r):
    ev = r["event_type"]

    if ev == "alt_acceptor":
        has_loss = pd.notna(r["diff_AL_ref_acceptor"]) and r["diff_AL_ref_acceptor"] == 0 and pd.notna(r["DS_AL"]) and r["DS_AL"] >= 0.8
        has_gain = pd.notna(r["diff_AG_obs_acceptor"]) and r["diff_AG_obs_acceptor"] == 0 and pd.notna(r["DS_AG"]) and r["DS_AG"] >= 0.5
        if has_loss and has_gain:
            return "very_strong_alt_acceptor_loss_plus_gain"
        elif has_loss:
            return "strong_alt_acceptor_ref_loss"
        elif has_gain:
            return "strong_alt_acceptor_observed_gain"

    if ev == "alt_donor":
        has_loss = pd.notna(r["diff_DL_ref_donor"]) and r["diff_DL_ref_donor"] == 0 and pd.notna(r["DS_DL"]) and r["DS_DL"] >= 0.8
        has_gain = pd.notna(r["diff_DG_obs_donor"]) and r["diff_DG_obs_donor"] == 0 and pd.notna(r["DS_DG"]) and r["DS_DG"] >= 0.5
        if has_loss and has_gain:
            return "very_strong_alt_donor_loss_plus_gain"
        elif has_loss:
            return "strong_alt_donor_ref_loss"
        elif has_gain:
            return "strong_alt_donor_observed_gain"

    if ev == "novel_exon":
        has_acc = pd.notna(r["diff_AG_obs_acceptor"]) and r["diff_AG_obs_acceptor"] == 0 and pd.notna(r["DS_AG"]) and r["DS_AG"] >= 0.5
        has_don = pd.notna(r["diff_DG_obs_donor"]) and r["diff_DG_obs_donor"] == 0 and pd.notna(r["DS_DG"]) and r["DS_DG"] >= 0.5
        if has_acc and has_don:
            return "very_strong_novel_exon_both_boundaries_supported"

    return "other"

strong["support_label"] = strong.apply(classify_row, axis=1)
very_strong["support_label"] = very_strong.apply(classify_row, axis=1)

#Best per PBID/event
best_key_cols = ["transcript_id", "event_type", "pb_exon_index", "observed_start", "observed_end"]

strong_best = (
    strong.sort_values(
        ["match_score", "maxDS", "DS_AL", "DS_DL", "DS_AG", "DS_DG"],
        ascending=[False, False, False, False, False, False]
    )
    .drop_duplicates(subset=best_key_cols, keep="first")
    .reset_index(drop=True)
)

#Summary by variant
summary_by_variant = (
    strong.groupby(
        ["variant_ID", "VariationID", "variant_Name", "VariantType"],
        dropna=False
    )
    .agg(
        n_strong_matches=("transcript_id", "size"),
        n_unique_pbids=("transcript_id", "nunique"),
        matched_pbids=("transcript_id", lambda x: ",".join(sorted(set(map(str, x))))),
        matched_event_types=("event_type", lambda x: ",".join(sorted(set(map(str, x))))),
        best_match_score=("match_score", "max"),
        best_maxDS=("maxDS", "max")
    )
    .reset_index()
    .sort_values(["n_unique_pbids", "best_match_score", "best_maxDS"], ascending=[False, False, False])
    .reset_index(drop=True)
)

#Reorder useful columns
preferred_cols = [
    "transcript_id", "cell_type", "event_type", "event_subtype",
    "pb_exon_index", "ref_exon_index",
    "observed_start", "observed_end",
    "observed_donor_pos", "observed_acceptor_pos",
    "reference_donor_pos", "reference_acceptor_pos",
    "variant_ID", "VariationID", "variant_Name", "VariantType",
    "DS_AG", "DS_AL", "DS_DG", "DS_DL",
    "diff_AG_obs_acceptor", "diff_DG_obs_donor",
    "diff_AL_ref_acceptor", "diff_DL_ref_donor",
    "score_AG_obs_acceptor", "score_DG_obs_donor",
    "score_AL_ref_acceptor", "score_DL_ref_donor",
    "match_score", "match_reason", "support_label"
]

strong = strong[[c for c in preferred_cols if c in strong.columns] + [c for c in strong.columns if c not in preferred_cols]]
very_strong = very_strong[[c for c in preferred_cols if c in very_strong.columns] + [c for c in very_strong.columns if c not in preferred_cols]]
strong_best = strong_best[[c for c in preferred_cols if c in strong_best.columns] + [c for c in strong_best.columns if c not in preferred_cols]]

#Save outputs
strong.to_csv(OUT_ALL, sep="\t", index=False)
very_strong.to_csv(OUT_VERY_STRONG, sep="\t", index=False)
strong_best.to_csv(OUT_BEST_PER_PBID, sep="\t", index=False)
summary_by_variant.to_csv(OUT_SUMMARY_BY_VARIANT, sep="\t", index=False)

with open(OUT_REPORT, "w") as out:
    out.write("ABCA4 strong SpliceAI support summary\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"Input file: {INFILE}\n\n")
    out.write(f"Total rows in input: {df.shape[0]}\n")
    out.write(f"Rows after excluding terminal exon events: {df_main.shape[0]}\n")
    out.write(f"Strong support rows: {strong.shape[0]}\n")
    out.write(f"Very strong support rows: {very_strong.shape[0]}\n")
    out.write(f"Best strong support per PBID/event: {strong_best.shape[0]}\n\n")

    out.write("Support label counts:\n")
    for k, v in strong["support_label"].value_counts().items():
        out.write(f"{k}\t{v}\n")
    out.write("\n")

    out.write("Top variants by number of unique PBIDs supported:\n")
    for _, r in summary_by_variant.head(20).iterrows():
        out.write(
            f"{r['variant_ID']}\t{r['VariationID']}\t{r['n_unique_pbids']}\t"
            f"{r['best_match_score']}\t{r['best_maxDS']}\t{r['variant_Name']}\n"
        )

print("Done.")
print("Saved:")
print(OUT_ALL)
print(OUT_VERY_STRONG)
print(OUT_BEST_PER_PBID)
print(OUT_SUMMARY_BY_VARIANT)
print(OUT_REPORT)

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################