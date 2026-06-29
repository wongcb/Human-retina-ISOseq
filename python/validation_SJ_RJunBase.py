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
rjun_file  = "detail_LS_annotation.txt"
out_fsm_ism_strict = "FSM_ISM_junction_validation_strict.csv"   # strict (0 bp, chr-clean)
out_nnc_nic_strict = "NNC_NIC_junction_validation_strict.csv"
out_fsm_ism_loose  = "FSM_ISM_junction_validation_loose.csv"    # loose (±1 bp, chr-clean)
out_nnc_nic_loose  = "NNC_NIC_junction_validation_loose.csv"
out_loose_summary  = "SJ_loose_match_summary.csv"               # summary for ±1/2/3 bp

#Read Iso-Seq classification
print("Reading Iso-Seq classification:", class_file)
df_class = pd.read_csv(class_file, sep="\t", na_values=["NA"])

#Keep only essential columns
df_class = df_class[["isoform", "FL", "structural_category"]].copy()
df_class["FL"] = pd.to_numeric(df_class["FL"], errors="coerce")

#Filter by FL >= 2 (FL2)
df_class = df_class[df_class["FL"] >= 2].copy()
print(f"Number of isoforms with FL >= 2: {df_class.shape[0]}")


def simplify_class(cat: str) -> str:
    """
    Map SQANTI3 structural_category to a simplified class label.
    """
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


#Map structural_category to FSM/ISM/NIC/NNC and keep only these four classes
df_class["sqanti_class"] = df_class["structural_category"].apply(simplify_class)
df_class = df_class[df_class["sqanti_class"].isin(["FSM", "ISM", "NIC", "NNC"])]

#Isoform ID sets for known-like vs novel classes
fsm_ism_ids = set(df_class[df_class["sqanti_class"].isin(["FSM", "ISM"])]["isoform"])
nnc_nic_ids = set(df_class[df_class["sqanti_class"].isin(["NIC", "NNC"])]["isoform"])

print(f"FSM/ISM isoforms: {len(fsm_ism_ids)}")
print(f"NNC/NIC isoforms: {len(nnc_nic_ids)}")

#Read Iso-Seq junction file
print("Reading Iso-Seq junction file:", junc_file)
df_junc = pd.read_csv(junc_file, sep="\t", na_values=["NA"])

#Basic columns from the junction file
df_junc["chrom_raw"]   = df_junc["chrom"].astype(str)
df_junc["start_iso"]   = df_junc["genomic_start_coord"].astype(int)
df_junc["end_iso"]     = df_junc["genomic_end_coord"].astype(int)
df_junc["strand_iso"]  = df_junc["strand"].astype(str)

#"chr1" -> "1", "chrX" -> "X"
df_junc["chrom_clean"] = df_junc["chrom_raw"].str.replace("^chr", "", regex=True)

#Subset junctions into FSM/ISM and NNC/NIC based on isoform IDs
junc_fsm_ism = df_junc[df_junc["isoform"].isin(fsm_ism_ids)].copy()
junc_nnc_nic = df_junc[df_junc["isoform"].isin(nnc_nic_ids)].copy()

print(f"FSM/ISM junctions: {len(junc_fsm_ism)}")
print(f"NNC/NIC junctions: {len(junc_nnc_nic)}")

#Read RJunBase and parse junction coordinates
print("Reading RJunBase file:", rjun_file)
rjun = pd.read_csv(rjun_file, sep="\t")

#Junction location: chromosome, start, end, and strand, and a cleaned chromosome.
loc = rjun["Junction location"].astype(str).str.extract(
    r'^(chr)?(?P<chrom>[^:]+):(?P<start>\d+)\|(?P<end>\d+):(?P<strand>[+-])',
    expand=True
)
rjun["chrom_clean"]  = loc["chrom"].astype(str)    # e.g., "1", "2", "X"
rjun["start_rjun"]   = loc["start"].astype(int)
rjun["end_rjun"]     = loc["end"].astype(int)
rjun["strand_rjun"]  = loc["strand"].astype(str)

#STRICT matching with chr-clean (0 bp tolerance)
print("Performing STRICT matching (0 bp, chr-clean)...")

#Build a chromosome-clean strict junction ID: "1:start|end:+"
rjun["sj_id_clean"] = (
    rjun["chrom_clean"].astype(str)
    + ":"
    + rjun["start_rjun"].astype(str)
    + "|"
    + rjun["end_rjun"].astype(str)
    + ":"
    + rjun["strand_rjun"].astype(str)
)
rjun_ids_strict = set(rjun["sj_id_clean"])

#Build the same style strict ID for Iso-Seq junctions (using chrom_clean)
junc_fsm_ism["sj_id_clean"] = (
    junc_fsm_ism["chrom_clean"].astype(str)
    + ":"
    + junc_fsm_ism["start_iso"].astype(str)
    + "|"
    + junc_fsm_ism["end_iso"].astype(str)
    + ":"
    + junc_fsm_ism["strand_iso"].astype(str)
)
junc_nnc_nic["sj_id_clean"] = (
    junc_nnc_nic["chrom_clean"].astype(str)
    + ":"
    + junc_nnc_nic["start_iso"].astype(str)
    + "|"
    + junc_nnc_nic["end_iso"].astype(str)
    + ":"
    + junc_nnc_nic["strand_iso"].astype(str)
)

#Strict in_RJunBase_strict (0 bp, chr-clean)
junc_fsm_ism["in_RJunBase_strict"] = junc_fsm_ism["sj_id_clean"].isin(rjun_ids_strict)
junc_nnc_nic["in_RJunBase_strict"] = junc_nnc_nic["sj_id_clean"].isin(rjun_ids_strict)

fsm_ism_strict_rate = junc_fsm_ism["in_RJunBase_strict"].mean()
nnc_nic_strict_rate = junc_nnc_nic["in_RJunBase_strict"].mean()

print(f"FSM/ISM strict match rate (0 bp, chr-clean): {fsm_ism_strict_rate:.3f}")
print(f"NNC/NIC strict match rate (0 bp, chr-clean): {nnc_nic_strict_rate:.3f}")

#LOOSE matching with ±1 / ±2 / ±3 bp tolerance
print("Preparing for LOOSE matching (±1 / ±2 / ±3 bp)...")

#Build an index for RJunBase for efficient coordinate-based matching:
#key = (chrom_clean, strand_rjun) -> array of [[start_rjun, end_rjun], ...]
rjun_index = {}
for (chrom, strand), sub in rjun.groupby(["chrom_clean", "strand_rjun"]):
    key = (str(chrom), str(strand))
    rjun_index[key] = sub[["start_rjun", "end_rjun"]].to_numpy()

print(f"Number of (chrom, strand) groups in RJunBase: {len(rjun_index)}")


def flag_loose_match(df: pd.DataFrame, tol: int = 1) -> pd.Series:
    """
    For each junction in df (FSM/ISM or NNC/NIC), check if there exists a
    junction in RJunBase on the same chromosome and strand whose start and end
    coordinates are within ±tol bp of the Iso-Seq junction.
    df must contain: 'chrom_clean', 'start_iso', 'end_iso', 'strand_iso'.
    """
    def _check(row):
        key = (str(row["chrom_clean"]), str(row["strand_iso"]))
        arr = rjun_index.get(key)
        if arr is None:
            # No junctions in RJunBase for this chromosome + strand
            return False
        s = int(row["start_iso"])
        e = int(row["end_iso"])
        # Vectorised comparison:
        # |start_rjun - s| <= tol AND |end_rjun - e| <= tol
        start_ok = np.abs(arr[:, 0] - s) <= tol
        end_ok   = np.abs(arr[:, 1] - e) <= tol
        return bool(np.any(start_ok & end_ok))

    return df.apply(_check, axis=1)


#Evaluate loose matching for multiple tolerances
tolerances = [1, 2, 3]
summary_rows = []

for tol in tolerances:
    print(f"\nPerforming LOOSE matching for FSM/ISM with tolerance ±{tol} bp...")
    col_fsm = f"in_RJunBase_loose_{tol}bp"
    junc_fsm_ism[col_fsm] = flag_loose_match(junc_fsm_ism, tol=tol)

    print(f"Performing LOOSE matching for NNC/NIC with tolerance ±{tol} bp...")
    col_nnc = f"in_RJunBase_loose_{tol}bp"
    junc_nnc_nic[col_nnc] = flag_loose_match(junc_nnc_nic, tol=tol)

    fsm_ism_loose_rate = junc_fsm_ism[col_fsm].mean()
    nnc_nic_loose_rate = junc_nnc_nic[col_nnc].mean()

    print(f"FSM/ISM loose match rate (±{tol} bp): {fsm_ism_loose_rate:.3f}")
    print(f"NNC/NIC loose match rate (±{tol} bp): {nnc_nic_loose_rate:.3f}")
    print("Number of matched junctions (loose) - FSM/ISM:",
          int(junc_fsm_ism[col_fsm].sum()))
    print("Number of matched junctions (loose) - NNC/NIC:",
          int(junc_nnc_nic[col_nnc].sum()))

    summary_rows.append({
        "tol_bp": tol,
        "FSM_ISM_rate": fsm_ism_loose_rate,
        "NNC_NIC_rate": nnc_nic_loose_rate,
        "FSM_ISM_matched": int(junc_fsm_ism[col_fsm].sum()),
        "NNC_NIC_matched": int(junc_nnc_nic[col_nnc].sum())
    })

#Save tolerance summary (±1/2/3 bp)
df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(out_loose_summary, index=False)
print("\n=== Summary of loose match rates (saved to {}) ===".format(out_loose_summary))
print(df_summary)

#For downstream convenience: use the ±1 bp result as the main "in_RJunBase_loose" column
junc_fsm_ism["in_RJunBase_loose"] = junc_fsm_ism["in_RJunBase_loose_1bp"]
junc_nnc_nic["in_RJunBase_loose"] = junc_nnc_nic["in_RJunBase_loose_1bp"]

#Export per-junction validation results
print("\nWriting strict match results (0 bp, chr-clean) to CSV...")
junc_fsm_ism[[
    "isoform", "chrom_clean", "start_iso", "end_iso",
    "strand_iso", "sj_id_clean", "in_RJunBase_strict"
]].to_csv(out_fsm_ism_strict, index=False)

junc_nnc_nic[[
    "isoform", "chrom_clean", "start_iso", "end_iso",
    "strand_iso", "sj_id_clean", "in_RJunBase_strict"
]].to_csv(out_nnc_nic_strict, index=False)

print("Writing loose match results (±1 bp, chr-clean) to CSV...")
junc_fsm_ism[[
    "isoform", "chrom_clean", "start_iso", "end_iso",
    "strand_iso", "in_RJunBase_loose"
]].to_csv(out_fsm_ism_loose, index=False)

junc_nnc_nic[[
    "isoform", "chrom_clean", "start_iso", "end_iso",
    "strand_iso", "in_RJunBase_loose"
]].to_csv(out_nnc_nic_loose, index=False)

print("Done.")

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################