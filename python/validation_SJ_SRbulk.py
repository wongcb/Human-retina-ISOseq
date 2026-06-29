########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np


#paths
class_file = "RetinaMerge_classification.filtered_lite_classification.txt"
junc_file  = "RetinaMerge_classification.filtered_lite_junctions.txt"
sr_sj_file = "merged_SJ_counts.tsv"
out_all_fsm_ism = "ALL_FSM_ISM_junction_shortread_support.csv"
out_all_nnc_nic = "ALL_NNC_NIC_junction_shortread_support.csv"
out_fl2_fsm_ism = "FL2_FSM_ISM_junction_shortread_support.csv"
out_fl2_nnc_nic = "FL2_NNC_NIC_junction_shortread_support.csv"
out_hc_fsm_ism = "HighConfident_FSM_ISM_junction_shortread_support.csv"
out_hc_nnc_nic = "HighConfident_NNC_NIC_junction_shortread_support.csv"
out_summary = "SJ_shortread_support_summary.csv"

#define functions
def simplify_class(cat: str) -> str:
    if cat == "full-splice_match":
        return "FSM"
    elif cat == "incomplete-splice_match":
        return "ISM"
    elif cat == "novel_in_catalog":
        return "NIC"
    elif cat == "novel_not_in_catalog":
        return "NNC"
    else:
        return "OTHER"


def is_true_series(x: pd.Series) -> pd.Series:
    """
    Convert possible TRUE/FALSE-like values into boolean TRUE flag.
    Handles TRUE, True, true, 1.
    """
    return (
        x.fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
        .isin(["TRUE", "T", "1"])
    )


def has_valid_polya_motif(x: pd.Series) -> pd.Series:
    """
    Valid polyA motif means not NA and not empty.
    This treats string values like 'NA', 'nan', 'None' as invalid.
    """
    s = (
        x.fillna("")
        .astype(str)
        .str.strip()
    )

    return (
        (s != "")
        & (~s.str.upper().isin(["NA", "NAN", "NONE", "NULL"]))
    )


def build_isoform_sets(df_class: pd.DataFrame):
    """
    Build isoform ID sets for ALL, FL2, and HighConfident,
    separated into FSM/ISM and NIC/NNC.

    HighConfident definition:
      FL >= 2
      within_cage_peak == TRUE
      polyA_motif is not NA and not empty
    """
    df_all = df_class.copy()

    df_fl2 = df_class[df_class["FL"] >= 2].copy()

    high_confident_mask = (
        (df_class["FL"] >= 2)
        & is_true_series(df_class["within_cage_peak"])
        & has_valid_polya_motif(df_class["polyA_motif"])
    )

    df_hc = df_class[high_confident_mask].copy()

    sets = {
        "ALL": {
            "FSM_ISM": set(df_all[df_all["sqanti_class"].isin(["FSM", "ISM"])]["isoform"]),
            "NNC_NIC": set(df_all[df_all["sqanti_class"].isin(["NIC", "NNC"])]["isoform"]),
        },
        "FL2": {
            "FSM_ISM": set(df_fl2[df_fl2["sqanti_class"].isin(["FSM", "ISM"])]["isoform"]),
            "NNC_NIC": set(df_fl2[df_fl2["sqanti_class"].isin(["NIC", "NNC"])]["isoform"]),
        },
        "HighConfident": {
            "FSM_ISM": set(df_hc[df_hc["sqanti_class"].isin(["FSM", "ISM"])]["isoform"]),
            "NNC_NIC": set(df_hc[df_hc["sqanti_class"].isin(["NIC", "NNC"])]["isoform"]),
        }
    }

    print("ALL classified isoforms:", df_all.shape[0])
    print("FL2 classified isoforms:", df_fl2.shape[0])
    print("HighConfident classified isoforms:", df_hc.shape[0])

    print("HighConfident FSM/ISM isoforms:", len(sets["HighConfident"]["FSM_ISM"]))
    print("HighConfident NNC/NIC isoforms:", len(sets["HighConfident"]["NNC_NIC"]))

    return sets


def build_sj_id(chrom, start, end, strand):
    return (
        chrom.astype(str)
        + ":"
        + start.astype(str)
        + "|"
        + end.astype(str)
        + ":"
        + strand.astype(str)
    )


def build_sr_exact_map(sr: pd.DataFrame) -> pd.DataFrame:
    """
    Build exact-match short-read support map keyed by sj_id_clean.
    """
    sr = sr.copy()

    sr["sj_id_clean"] = build_sj_id(
        sr["chrom_clean"], sr["start"], sr["end"], sr["strand"]
    )

    sr_map = sr[[
        "sj_id_clean",
        "unique_reads_sum",
        "multi_reads_sum",
        "samples_detected",
        "intron_motif",
        "annotated"
    ]].drop_duplicates().copy()

    sr_map = sr_map.rename(columns={
        "intron_motif": "intron_motif_sr",
        "annotated": "annotated_sr"
    })

    return sr_map


def build_sr_loose_index(sr: pd.DataFrame):
    """
    Build index for loose matching:
    key = (chrom_clean, strand)
    value = ndarray of [start, end, unique_reads_sum, multi_reads_sum, samples_detected]
    """
    loose_index = {}

    for (chrom, strand), sub in sr.groupby(["chrom_clean", "strand"], sort=False):
        key = (str(chrom), str(strand))
        loose_index[key] = sub[[
            "start",
            "end",
            "unique_reads_sum",
            "multi_reads_sum",
            "samples_detected"
        ]].to_numpy()

    return loose_index


def annotate_strict_support(df_junc_sub: pd.DataFrame, sr_map: pd.DataFrame) -> pd.DataFrame:
    """
    Strict exact match on chrom_clean, start, end, strand.
    """
    out = df_junc_sub.copy()

    out["sj_id_clean"] = build_sj_id(
        out["chrom_clean"],
        out["start_iso"],
        out["end_iso"],
        out["strand_iso"]
    )

    out = out.merge(sr_map, how="left", on="sj_id_clean")

    out["in_shortread_strict"] = ~out["unique_reads_sum"].isna()

    fill_cols = [
        "unique_reads_sum",
        "multi_reads_sum",
        "samples_detected",
        "intron_motif_sr",
        "annotated_sr"
    ]

    for c in fill_cols:
        if c in out.columns:
            out[c] = out[c].fillna(0)

    return out


def annotate_loose_support(df: pd.DataFrame, loose_index: dict, tol: int = 1) -> pd.DataFrame:
    """
    Add ±tol bp loose support columns.

    Matching rule:
      same chrom_clean,
      same strand,
      abs(start_sr - start_iso) <= tol,
      abs(end_sr - end_iso) <= tol

    For loose support evidence columns:
      - max unique_reads_sum among matched short-read junctions
      - max multi_reads_sum among matched short-read junctions
      - max samples_detected among matched short-read junctions
    """
    out = df.copy()

    loose_flag = []
    loose_unique = []
    loose_multi = []
    loose_samples = []

    for _, row in out.iterrows():
        key = (str(row["chrom_clean"]), str(row["strand_iso"]))
        arr = loose_index.get(key)

        if arr is None or len(arr) == 0:
            loose_flag.append(False)
            loose_unique.append(0)
            loose_multi.append(0)
            loose_samples.append(0)
            continue

        s = int(row["start_iso"])
        e = int(row["end_iso"])

        start_ok = np.abs(arr[:, 0].astype(int) - s) <= tol
        end_ok = np.abs(arr[:, 1].astype(int) - e) <= tol
        mask = start_ok & end_ok

        if np.any(mask):
            matched = arr[mask]
            loose_flag.append(True)
            loose_unique.append(int(np.max(matched[:, 2].astype(int))))
            loose_multi.append(int(np.max(matched[:, 3].astype(int))))
            loose_samples.append(int(np.max(matched[:, 4].astype(int))))
        else:
            loose_flag.append(False)
            loose_unique.append(0)
            loose_multi.append(0)
            loose_samples.append(0)

    out[f"in_shortread_loose_{tol}bp"] = loose_flag
    out[f"loose_{tol}bp_unique_reads_max"] = loose_unique
    out[f"loose_{tol}bp_multi_reads_max"] = loose_multi
    out[f"loose_{tol}bp_samples_detected_max"] = loose_samples

    return out


def summarise_group(df: pd.DataFrame, set_name: str, group_name: str) -> dict:
    """
    Build summary row for one group.
    """
    n = len(df)

    n_supported_strict = int(df["in_shortread_strict"].sum()) if n > 0 else 0
    strict_rate = n_supported_strict / n if n > 0 else np.nan

    n_supported_loose = int(df["in_shortread_loose_1bp"].sum()) if n > 0 else 0
    loose_rate = n_supported_loose / n if n > 0 else np.nan

    strict_supported = df[df["in_shortread_strict"]].copy() if n > 0 else df.copy()
    loose_supported = df[df["in_shortread_loose_1bp"]].copy() if n > 0 else df.copy()

    return {
        "set_name": set_name,
        "group": group_name,
        "n_junctions": n,

        "n_supported_strict": n_supported_strict,
        "support_rate_strict": strict_rate,

        "n_supported_loose_1bp": n_supported_loose,
        "support_rate_loose_1bp": loose_rate,

        "median_unique_reads_sum_strict_supported":
            strict_supported["unique_reads_sum"].median() if len(strict_supported) > 0 else np.nan,

        "median_samples_detected_strict_supported":
            strict_supported["samples_detected"].median() if len(strict_supported) > 0 else np.nan,

        "median_loose_1bp_unique_reads_max_supported":
            loose_supported["loose_1bp_unique_reads_max"].median() if len(loose_supported) > 0 else np.nan,

        "median_loose_1bp_samples_detected_max_supported":
            loose_supported["loose_1bp_samples_detected_max"].median() if len(loose_supported) > 0 else np.nan,

        "n_supported_unique_ge_1_strict":
            int((df["unique_reads_sum"] >= 1).sum()) if n > 0 else 0,

        "n_supported_unique_ge_3_strict":
            int((df["unique_reads_sum"] >= 3).sum()) if n > 0 else 0,

        "n_supported_unique_ge_5_strict":
            int((df["unique_reads_sum"] >= 5).sum()) if n > 0 else 0,

        "n_supported_samples_ge_2_strict":
            int((df["samples_detected"] >= 2).sum()) if n > 0 else 0,

        "n_supported_loose_unique_ge_1":
            int((df["loose_1bp_unique_reads_max"] >= 1).sum()) if n > 0 else 0,

        "n_supported_loose_unique_ge_3":
            int((df["loose_1bp_unique_reads_max"] >= 3).sum()) if n > 0 else 0,

        "n_supported_loose_samples_ge_2":
            int((df["loose_1bp_samples_detected_max"] >= 2).sum()) if n > 0 else 0,
    }


def process_one_group(
    df_junc: pd.DataFrame,
    isoform_set: set,
    sr_map: pd.DataFrame,
    loose_index: dict,
    set_name: str,
    group_name: str
):
    """
    Subset junctions by isoform set, annotate strict/loose support,
    and return per-junction result plus summary row.
    """
    print(f"\nProcessing {set_name} {group_name}...")

    df_sub = df_junc[df_junc["isoform"].isin(isoform_set)].copy()

    df_sub = annotate_strict_support(df_sub, sr_map)
    df_sub = annotate_loose_support(df_sub, loose_index, tol=1)

    summary = summarise_group(df_sub, set_name, group_name)

    print(f"{set_name} {group_name} junctions:", summary["n_junctions"])
    print(f"{set_name} {group_name} strict support rate:", summary["support_rate_strict"])
    print(f"{set_name} {group_name} loose ±1bp support rate:", summary["support_rate_loose_1bp"])

    return df_sub, summary


#Read classification
print("Reading classification:", class_file)
df_class = pd.read_csv(class_file, sep="\t", na_values=["NA"], low_memory=False)

required_class_cols = [
    "isoform",
    "FL",
    "structural_category",
    "within_cage_peak",
    "polyA_motif"
]

missing_class_cols = [c for c in required_class_cols if c not in df_class.columns]
if len(missing_class_cols) > 0:
    raise ValueError(
        "Missing required columns in classification file: "
        + ", ".join(missing_class_cols)
    )

df_class = df_class[required_class_cols].copy()

df_class["FL"] = pd.to_numeric(df_class["FL"], errors="coerce")
df_class["sqanti_class"] = df_class["structural_category"].apply(simplify_class)

df_class = df_class[df_class["sqanti_class"].isin(["FSM", "ISM", "NIC", "NNC"])].copy()

print("Total classified isoforms (FSM/ISM/NIC/NNC):", df_class.shape[0])
print("FL>=2 classified isoforms (FSM/ISM/NIC/NNC):", int((df_class["FL"] >= 2).sum()))

hc_mask_check = (
    (df_class["FL"] >= 2)
    & is_true_series(df_class["within_cage_peak"])
    & has_valid_polya_motif(df_class["polyA_motif"])
)

print(
    "HighConfident classified isoforms before FSM/ISM/NIC/NNC grouping:",
    int(hc_mask_check.sum())
)

isoform_sets = build_isoform_sets(df_class)


#Read Iso-Seq junction file
print("\nReading Iso-Seq junction file:", junc_file)
df_junc = pd.read_csv(junc_file, sep="\t", na_values=["NA"], low_memory=False)

required_junc_cols = [
    "isoform",
    "chrom",
    "strand",
    "genomic_start_coord",
    "genomic_end_coord"
]

missing_junc_cols = [c for c in required_junc_cols if c not in df_junc.columns]
if len(missing_junc_cols) > 0:
    raise ValueError(
        "Missing required columns in junction file: "
        + ", ".join(missing_junc_cols)
    )

df_junc["chrom_raw"] = df_junc["chrom"].astype(str)
df_junc["chrom_clean"] = df_junc["chrom_raw"].str.replace("^chr", "", regex=True)
df_junc["start_iso"] = pd.to_numeric(df_junc["genomic_start_coord"], errors="raise").astype(int)
df_junc["end_iso"] = pd.to_numeric(df_junc["genomic_end_coord"], errors="raise").astype(int)
df_junc["strand_iso"] = df_junc["strand"].astype(str)

print("Total Iso-Seq junction rows:", df_junc.shape[0])


#Read merged short-read SJ table
print("\nReading merged short-read SJ table:", sr_sj_file)
sr = pd.read_csv(sr_sj_file, sep="\t", na_values=["NA"], low_memory=False)

required_sr_cols = [
    "chr",
    "start",
    "end",
    "strand",
    "intron_motif",
    "annotated",
    "unique_reads_sum",
    "multi_reads_sum",
    "samples_detected"
]

missing_sr_cols = [c for c in required_sr_cols if c not in sr.columns]
if len(missing_sr_cols) > 0:
    raise ValueError(
        "Missing required columns in short-read SJ file: "
        + ", ".join(missing_sr_cols)
    )

sr = sr[~((sr["chr"] == "chr") & (sr["start"].astype(str) == "start"))].copy()

sr["chr"] = sr["chr"].astype(str)
sr["chrom_clean"] = sr["chr"].str.replace("^chr", "", regex=True)

sr["start"] = pd.to_numeric(sr["start"], errors="raise").astype(int)
sr["end"] = pd.to_numeric(sr["end"], errors="raise").astype(int)

#STAR strand encoding: 1=+, 2=-, 0=undefined
sr["strand"] = sr["strand"].astype(str)
sr.loc[sr["strand"] == "1", "strand"] = "+"
sr.loc[sr["strand"] == "2", "strand"] = "-"
sr.loc[sr["strand"] == "0", "strand"] = "."

sr["intron_motif"] = pd.to_numeric(sr["intron_motif"], errors="raise").astype(int)
sr["annotated"] = pd.to_numeric(sr["annotated"], errors="raise").astype(int)
sr["unique_reads_sum"] = pd.to_numeric(sr["unique_reads_sum"], errors="raise").astype(int)
sr["multi_reads_sum"] = pd.to_numeric(sr["multi_reads_sum"], errors="raise").astype(int)
sr["samples_detected"] = pd.to_numeric(sr["samples_detected"], errors="raise").astype(int)

sr_map = build_sr_exact_map(sr)
sr_ids = set(sr_map["sj_id_clean"])
loose_index = build_sr_loose_index(sr)

print("Number of merged short-read unique junctions:", len(sr_ids))


#Process ALL / FL2 / HighConfident
summary_rows = []

print("\n=== Processing ALL isoforms ===")
all_fsm_ism, summary = process_one_group(
    df_junc=df_junc,
    isoform_set=isoform_sets["ALL"]["FSM_ISM"],
    sr_map=sr_map,
    loose_index=loose_index,
    set_name="ALL",
    group_name="FSM_ISM"
)
summary_rows.append(summary)

all_nnc_nic, summary = process_one_group(
    df_junc=df_junc,
    isoform_set=isoform_sets["ALL"]["NNC_NIC"],
    sr_map=sr_map,
    loose_index=loose_index,
    set_name="ALL",
    group_name="NNC_NIC"
)
summary_rows.append(summary)


print("\n=== Processing FL2 isoforms ===")
fl2_fsm_ism, summary = process_one_group(
    df_junc=df_junc,
    isoform_set=isoform_sets["FL2"]["FSM_ISM"],
    sr_map=sr_map,
    loose_index=loose_index,
    set_name="FL2",
    group_name="FSM_ISM"
)
summary_rows.append(summary)

fl2_nnc_nic, summary = process_one_group(
    df_junc=df_junc,
    isoform_set=isoform_sets["FL2"]["NNC_NIC"],
    sr_map=sr_map,
    loose_index=loose_index,
    set_name="FL2",
    group_name="NNC_NIC"
)
summary_rows.append(summary)


print("\n=== Processing HighConfident isoforms ===")
print("Definition: FL >= 2, within_cage_peak == TRUE, polyA_motif not empty/NA")

hc_fsm_ism, summary = process_one_group(
    df_junc=df_junc,
    isoform_set=isoform_sets["HighConfident"]["FSM_ISM"],
    sr_map=sr_map,
    loose_index=loose_index,
    set_name="HighConfident",
    group_name="FSM_ISM"
)
summary_rows.append(summary)

hc_nnc_nic, summary = process_one_group(
    df_junc=df_junc,
    isoform_set=isoform_sets["HighConfident"]["NNC_NIC"],
    sr_map=sr_map,
    loose_index=loose_index,
    set_name="HighConfident",
    group_name="NNC_NIC"
)
summary_rows.append(summary)


#Save outputs
print("\nWriting per-junction output files...")

all_fsm_ism.to_csv(out_all_fsm_ism, index=False)
all_nnc_nic.to_csv(out_all_nnc_nic, index=False)

fl2_fsm_ism.to_csv(out_fl2_fsm_ism, index=False)
fl2_nnc_nic.to_csv(out_fl2_nnc_nic, index=False)

hc_fsm_ism.to_csv(out_hc_fsm_ism, index=False)
hc_nnc_nic.to_csv(out_hc_nnc_nic, index=False)

df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(out_summary, index=False)

print("\n=== Summary ===")
print(df_summary)

print("\nOutput files:")
print(out_all_fsm_ism)
print(out_all_nnc_nic)
print(out_fl2_fsm_ism)
print(out_fl2_nnc_nic)
print(out_hc_fsm_ism)
print(out_hc_nnc_nic)
print(out_summary)

print("\nDone.")

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################