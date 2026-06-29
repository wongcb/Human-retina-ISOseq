########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import math
import pandas as pd

#paths
EVENTS_TSV = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/ABCA4_spliceai_ready/abca4_spliceai_comparison_events.tsv"
SPLICEAI_VCF = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/clinvar_spliceai_output/clinvar_ABCA4.spliceai_output.vcf"

OUT_DIR = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/ABCA4_spliceai_match"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_ALL_MATCHES = os.path.join(OUT_DIR, "abca4_spliceai_variant_event_matches.tsv")
OUT_BEST_MATCHES = os.path.join(OUT_DIR, "abca4_spliceai_best_matches.tsv")
OUT_VARIANT_SUMMARY = os.path.join(OUT_DIR, "abca4_spliceai_variant_summary.tsv")
OUT_REPORT = os.path.join(OUT_DIR, "abca4_spliceai_match_report.txt")


#parameters
MAX_POS_DIFF_BP = 25
MIN_SCORE_WEAK = 0.05
MIN_SCORE_MODERATE = 0.20
MIN_SCORE_STRONG = 0.50


#define functions
def safe_int(x):
    if pd.isna(x):
        return None
    try:
        if str(x).strip() == "":
            return None
        return int(float(x))
    except Exception:
        return None


def safe_float(x):
    if pd.isna(x):
        return None
    try:
        if str(x).strip() == "":
            return None
        return float(x)
    except Exception:
        return None


def absdiff(a, b):
    if a is None or b is None:
        return None
    return abs(int(a) - int(b))


def parse_info_field(info_str):
    d = {}
    if pd.isna(info_str):
        return d
    for item in str(info_str).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def parse_spliceai_annotation(spliceai_str):
    """
    Format:
    ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

    SpliceAI field can contain multiple comma-separated annotations,
    but for this project we keep entries where SYMBOL == ABCA4 and
    choose the one with the largest max delta score.
    """
    anns = []
    if spliceai_str is None:
        return anns

    for token in str(spliceai_str).split(","):
        parts = token.split("|")
        if len(parts) != 10:
            continue

        ann = {
            "ALLELE": parts[0],
            "SYMBOL": parts[1],
            "DS_AG": safe_float(parts[2]),
            "DS_AL": safe_float(parts[3]),
            "DS_DG": safe_float(parts[4]),
            "DS_DL": safe_float(parts[5]),
            "DP_AG": safe_int(parts[6]),
            "DP_AL": safe_int(parts[7]),
            "DP_DG": safe_int(parts[8]),
            "DP_DL": safe_int(parts[9]),
        }
        ann["maxDS"] = max([
            x for x in [ann["DS_AG"], ann["DS_AL"], ann["DS_DG"], ann["DS_DL"]] if x is not None
        ] + [0.0])
        anns.append(ann)

    return anns


def choose_best_abca4_annotation(anns):
    anns2 = [a for a in anns if a.get("SYMBOL") == "ABCA4"]
    if not anns2:
        return None
    anns2 = sorted(anns2, key=lambda x: x["maxDS"], reverse=True)
    return anns2[0]


def score_to_label(x):
    if x is None:
        return "NA"
    if x >= MIN_SCORE_STRONG:
        return "strong"
    if x >= MIN_SCORE_MODERATE:
        return "moderate"
    if x >= MIN_SCORE_WEAK:
        return "weak"
    return "very_weak"


def compute_predicted_positions(pos, ann):
    """
    SpliceAI predicted genomic positions = variant POS + DP_*
    """
    out = {}
    for k in ["AG", "AL", "DG", "DL"]:
        dp = ann.get(f"DP_{k}")
        out[f"pred_{k}_pos"] = None if dp is None else int(pos) + int(dp)
    return out


def closeness_score(diff_bp, max_bp=MAX_POS_DIFF_BP):
    if diff_bp is None:
        return 0.0
    if diff_bp > max_bp:
        return 0.0
    return 1.0 - (diff_bp / max_bp)


def combined_score(delta_score, diff_bp, max_bp=MAX_POS_DIFF_BP):
    """
    Conservative combined score:
    - if no positional agreement, no match
    - otherwise weight both SpliceAI strength and positional closeness
    """
    if delta_score is None:
        return 0.0
    c = closeness_score(diff_bp, max_bp=max_bp)
    if c <= 0:
        return 0.0
    return float(delta_score) * c


#Load events
events = pd.read_csv(EVENTS_TSV, sep="\t")
for col in [
    "pb_exon_index", "ref_exon_index",
    "observed_start", "observed_end",
    "observed_donor_pos", "observed_acceptor_pos",
    "reference_donor_pos", "reference_acceptor_pos"
]:
    if col in events.columns:
        events[col] = events[col].apply(safe_int)

if events.empty:
    raise ValueError("Event table is empty.")


#Parse SpliceAI VCF
records = []

with open(SPLICEAI_VCF, "r") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
        pos = int(pos)
        infod = parse_info_field(info)
        anns = parse_spliceai_annotation(infod.get("SpliceAI"))
        ann = choose_best_abca4_annotation(anns)

        if ann is None:
            continue

        pred = compute_predicted_positions(pos, ann)

        records.append({
            "CHROM": chrom,
            "POS": pos,
            "ID": vid,
            "REF": ref,
            "ALT": alt,
            "VariationID": infod.get("VariationID"),
            "Name": infod.get("Name"),
            "VariantType": infod.get("VariantType"),
            "DS_AG": ann["DS_AG"],
            "DS_AL": ann["DS_AL"],
            "DS_DG": ann["DS_DG"],
            "DS_DL": ann["DS_DL"],
            "DP_AG": ann["DP_AG"],
            "DP_AL": ann["DP_AL"],
            "DP_DG": ann["DP_DG"],
            "DP_DL": ann["DP_DL"],
            "pred_AG_pos": pred["pred_AG_pos"],
            "pred_AL_pos": pred["pred_AL_pos"],
            "pred_DG_pos": pred["pred_DG_pos"],
            "pred_DL_pos": pred["pred_DL_pos"],
            "maxDS": ann["maxDS"],
        })

spliceai_df = pd.DataFrame(records)
if spliceai_df.empty:
    raise ValueError("No ABCA4 SpliceAI annotations parsed from VCF.")

#deduplicate exact repeated entries
spliceai_df = spliceai_df.drop_duplicates(
    subset=["CHROM", "POS", "REF", "ALT", "ID", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]
).reset_index(drop=True)


#Matching
match_rows = []

for _, ev in events.iterrows():
    ev_type = str(ev["event_type"])
    ev_chr = ev["seqnames"]
    ev_obs_donor = safe_int(ev.get("observed_donor_pos"))
    ev_obs_acceptor = safe_int(ev.get("observed_acceptor_pos"))
    ev_ref_donor = safe_int(ev.get("reference_donor_pos"))
    ev_ref_acceptor = safe_int(ev.get("reference_acceptor_pos"))

    #restrict to same chromosome
    cand = spliceai_df.loc[spliceai_df["CHROM"] == ev_chr].copy()
    if cand.empty:
        continue

    for _, var in cand.iterrows():
        reason = []
        match_score = 0.0

        #precompute diffs
        diff_AG_obsA = absdiff(var["pred_AG_pos"], ev_obs_acceptor)
        diff_DG_obsD = absdiff(var["pred_DG_pos"], ev_obs_donor)
        diff_AL_refA = absdiff(var["pred_AL_pos"], ev_ref_acceptor)
        diff_DL_refD = absdiff(var["pred_DL_pos"], ev_ref_donor)

        diff_AG_refA = absdiff(var["pred_AG_pos"], ev_ref_acceptor)
        diff_DG_refD = absdiff(var["pred_DG_pos"], ev_ref_donor)
        diff_AL_obsA = absdiff(var["pred_AL_pos"], ev_obs_acceptor)
        diff_DL_obsD = absdiff(var["pred_DL_pos"], ev_obs_donor)

        score_AG_obsA = combined_score(var["DS_AG"], diff_AG_obsA)
        score_DG_obsD = combined_score(var["DS_DG"], diff_DG_obsD)
        score_AL_refA = combined_score(var["DS_AL"], diff_AL_refA)
        score_DL_refD = combined_score(var["DS_DL"], diff_DL_refD)

        score_AG_refA = combined_score(var["DS_AG"], diff_AG_refA)
        score_DG_refD = combined_score(var["DS_DG"], diff_DG_refD)
        score_AL_obsA = combined_score(var["DS_AL"], diff_AL_obsA)
        score_DL_obsD = combined_score(var["DS_DL"], diff_DL_obsD)

        if ev_type == "novel_exon":
            #best support = gained acceptor near observed acceptor
            #plus gained donor near observed donor
            match_score = (score_AG_obsA + score_DG_obsD) / 2.0
            if score_AG_obsA > 0:
                reason.append(f"AG->obs_acceptor(diff={diff_AG_obsA},DS={var['DS_AG']})")
            if score_DG_obsD > 0:
                reason.append(f"DG->obs_donor(diff={diff_DG_obsD},DS={var['DS_DG']})")

        elif ev_type == "alt_acceptor":
            #either gain at observed acceptor or loss at reference acceptor
            match_score = max(score_AG_obsA, score_AL_refA, score_AG_refA, score_AL_obsA)
            if score_AG_obsA == match_score and score_AG_obsA > 0:
                reason.append(f"AG->obs_acceptor(diff={diff_AG_obsA},DS={var['DS_AG']})")
            if score_AL_refA == match_score and score_AL_refA > 0:
                reason.append(f"AL->ref_acceptor(diff={diff_AL_refA},DS={var['DS_AL']})")
            if score_AG_refA == match_score and score_AG_refA > 0:
                reason.append(f"AG->ref_acceptor(diff={diff_AG_refA},DS={var['DS_AG']})")
            if score_AL_obsA == match_score and score_AL_obsA > 0:
                reason.append(f"AL->obs_acceptor(diff={diff_AL_obsA},DS={var['DS_AL']})")

        elif ev_type == "alt_donor":
            #either gain at observed donor or loss at reference donor
            match_score = max(score_DG_obsD, score_DL_refD, score_DG_refD, score_DL_obsD)
            if score_DG_obsD == match_score and score_DG_obsD > 0:
                reason.append(f"DG->obs_donor(diff={diff_DG_obsD},DS={var['DS_DG']})")
            if score_DL_refD == match_score and score_DL_refD > 0:
                reason.append(f"DL->ref_donor(diff={diff_DL_refD},DS={var['DS_DL']})")
            if score_DG_refD == match_score and score_DG_refD > 0:
                reason.append(f"DG->ref_donor(diff={diff_DG_refD},DS={var['DS_DG']})")
            if score_DL_obsD == match_score and score_DL_obsD > 0:
                reason.append(f"DL->obs_donor(diff={diff_DL_obsD},DS={var['DS_DL']})")

        elif ev_type == "alternative_first_exon":
            #weaker heuristic
            match_score = max(score_AG_obsA, score_DG_obsD, score_AL_refA, score_DL_refD)
            if match_score > 0:
                reason.append("terminal_exon_boundary_match")

        elif ev_type == "alternative_last_exon":
            match_score = max(score_AG_obsA, score_DG_obsD, score_AL_refA, score_DL_refD)
            if match_score > 0:
                reason.append("terminal_exon_boundary_match")

        elif ev_type == "partial_intron_retention":
            #nearby loss/gain around retained region boundary
            match_score = max(score_AG_obsA, score_DG_obsD, score_AL_refA, score_DL_refD)
            if match_score > 0:
                reason.append("partial_intron_retention_boundary_match")

        elif ev_type == "intron_retention":
            match_score = max(score_AL_refA, score_DL_refD, score_AG_obsA, score_DG_obsD)
            if match_score > 0:
                reason.append("intron_retention_boundary_match")

        elif ev_type == "exon_skipping":
            match_score = max(score_AL_refA, score_DL_refD, score_AG_obsA, score_DG_obsD)
            if match_score > 0:
                reason.append("exon_skipping_boundary_match")

        else:
            #generic fallback
            match_score = max(score_AG_obsA, score_DG_obsD, score_AL_refA, score_DL_refD)

        if match_score <= 0:
            continue

        match_rows.append({
            "transcript_id": ev["transcript_id"],
            "cell_type": ev.get("cell_type"),
            "event_type": ev_type,
            "event_subtype": ev.get("event_subtype"),
            "pb_exon_index": ev.get("pb_exon_index"),
            "ref_exon_index": ev.get("ref_exon_index"),
            "event_seqnames": ev["seqnames"],
            "event_strand": ev.get("strand"),
            "observed_start": ev.get("observed_start"),
            "observed_end": ev.get("observed_end"),
            "observed_donor_pos": ev_obs_donor,
            "observed_acceptor_pos": ev_obs_acceptor,
            "reference_donor_pos": ev_ref_donor,
            "reference_acceptor_pos": ev_ref_acceptor,
            "affected_region": ev.get("affected_region"),
            "spliceai_expected_signal": ev.get("spliceai_expected_signal"),

            "variant_CHROM": var["CHROM"],
            "variant_POS": var["POS"],
            "variant_ID": var["ID"],
            "VariationID": var["VariationID"],
            "variant_Name": var["Name"],
            "variant_REF": var["REF"],
            "variant_ALT": var["ALT"],
            "VariantType": var["VariantType"],

            "DS_AG": var["DS_AG"],
            "DS_AL": var["DS_AL"],
            "DS_DG": var["DS_DG"],
            "DS_DL": var["DS_DL"],
            "DP_AG": var["DP_AG"],
            "DP_AL": var["DP_AL"],
            "DP_DG": var["DP_DG"],
            "DP_DL": var["DP_DL"],
            "pred_AG_pos": var["pred_AG_pos"],
            "pred_AL_pos": var["pred_AL_pos"],
            "pred_DG_pos": var["pred_DG_pos"],
            "pred_DL_pos": var["pred_DL_pos"],
            "maxDS": var["maxDS"],
            "maxDS_label": score_to_label(var["maxDS"]),

            "diff_AG_obs_acceptor": diff_AG_obsA,
            "diff_DG_obs_donor": diff_DG_obsD,
            "diff_AL_ref_acceptor": diff_AL_refA,
            "diff_DL_ref_donor": diff_DL_refD,

            "score_AG_obs_acceptor": score_AG_obsA,
            "score_DG_obs_donor": score_DG_obsD,
            "score_AL_ref_acceptor": score_AL_refA,
            "score_DL_ref_donor": score_DL_refD,

            "match_score": match_score,
            "match_reason": "; ".join(reason) if reason else "generic"
        })

all_matches = pd.DataFrame(match_rows)

if all_matches.empty:
    raise ValueError(
        "No matches found. This can happen if all SpliceAI predicted positions are too far from event coordinates "
        f"given MAX_POS_DIFF_BP={MAX_POS_DIFF_BP}. Try increasing it to 50 or 100."
    )

all_matches = all_matches.sort_values(
    ["match_score", "maxDS", "transcript_id", "event_type"],
    ascending=[False, False, True, True]
).reset_index(drop=True)


#Best match per isoform
event_key_cols = [
    "transcript_id", "event_type", "event_subtype",
    "pb_exon_index", "observed_start", "observed_end"
]

best_matches = (
    all_matches
    .sort_values(["match_score", "maxDS"], ascending=[False, False])
    .drop_duplicates(subset=event_key_cols, keep="first")
    .reset_index(drop=True)
)

#Variant summary
variant_summary = (
    all_matches
    .groupby(["variant_ID", "VariationID", "variant_Name", "VariantType"], dropna=False)
    .agg(
        n_matches=("match_score", "size"),
        best_match_score=("match_score", "max"),
        best_maxDS=("maxDS", "max"),
        matched_transcripts=("transcript_id", lambda x: ",".join(sorted(set(map(str, x))))),
        matched_event_types=("event_type", lambda x: ",".join(sorted(set(map(str, x)))))
    )
    .reset_index()
    .sort_values(["best_match_score", "best_maxDS"], ascending=[False, False])
    .reset_index(drop=True)
)

#Save outputs
all_matches.to_csv(OUT_ALL_MATCHES, sep="\t", index=False)
best_matches.to_csv(OUT_BEST_MATCHES, sep="\t", index=False)
variant_summary.to_csv(OUT_VARIANT_SUMMARY, sep="\t", index=False)

with open(OUT_REPORT, "w") as out:
    out.write("ABCA4 SpliceAI vs Iso-Seq event matching report\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"EVENTS_TSV = {EVENTS_TSV}\n")
    out.write(f"SPLICEAI_VCF = {SPLICEAI_VCF}\n\n")
    out.write(f"Events loaded: {events.shape[0]}\n")
    out.write(f"SpliceAI variants parsed: {spliceai_df.shape[0]}\n")
    out.write(f"All matches: {all_matches.shape[0]}\n")
    out.write(f"Best matches: {best_matches.shape[0]}\n\n")

    out.write("Event type counts in best matches:\n")
    for k, v in best_matches["event_type"].value_counts().items():
        out.write(f"{k}\t{v}\n")
    out.write("\n")

    out.write("Top 20 variants by best match score:\n")
    for _, r in variant_summary.head(20).iterrows():
        out.write(
            f"{r['variant_ID']}\t{r['VariationID']}\t{r['best_match_score']:.4f}\t"
            f"{r['best_maxDS']:.4f}\t{r['matched_event_types']}\n"
        )

print("Done.")
print("Saved:")
print(OUT_ALL_MATCHES)
print(OUT_BEST_MATCHES)
print(OUT_VARIANT_SUMMARY)
print(OUT_REPORT)

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################