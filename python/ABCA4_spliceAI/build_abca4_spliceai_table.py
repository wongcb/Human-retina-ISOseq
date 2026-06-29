########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import gzip
from collections import Counter

import pandas as pd

#paths
REF_GTF = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/GRCh38p14/GRCh38.p14.v44/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz"

ABCA4_EXPORT_DIR = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/ABCA4_export"
EXTRACTED_DIR = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/extracted_pbids"

EXON_STATUS_CSV = os.path.join(ABCA4_EXPORT_DIR, "res_ABCA4_exon_status.csv")
JUNCTION_CSV = os.path.join(ABCA4_EXPORT_DIR, "ABCA4_all_junctions.csv")
GFF_FILE = os.path.join(EXTRACTED_DIR, "selected_pbids.gff")
FASTA_FILE = os.path.join(EXTRACTED_DIR, "selected_pbids.fa")

OUT_DIR = os.path.join("/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI", "ABCA4_spliceai_ready")
os.makedirs(OUT_DIR, exist_ok=True)

#define functions
def coord_key(start, end):
    return f"{int(start)}-{int(end)}"


def parse_gtf_attributes(attr_str):
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


def parse_gff_attributes(attr_str):
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


def safe_mode(vals):
    vals = [x for x in vals if pd.notna(x)]
    if not vals:
        return None
    return Counter(vals).most_common(1)[0][0]


def ordered_exons(df, strand):
    d = df.sort_values(["start", "end"]).copy()
    if strand == "-":
        d = d.iloc[::-1].copy()
    return d.reset_index(drop=True)


def exon_overlap_len(a_start, a_end, b_start, b_end):
    x1 = max(int(a_start), int(b_start))
    x2 = min(int(a_end), int(b_end))
    return max(0, x2 - x1 + 1)


def intervals_overlap(a_start, a_end, b_start, b_end):
    return exon_overlap_len(a_start, a_end, b_start, b_end) > 0


def parse_fasta_lengths(fasta_file):
    seq_lens = {}
    curr_id = None
    curr_len = 0
    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_id is not None:
                    seq_lens[curr_id] = curr_len
                curr_id = line[1:].split()[0]
                curr_len = 0
            else:
                curr_len += len(line)
        if curr_id is not None:
            seq_lens[curr_id] = curr_len
    return seq_lens


#Load inputs
for fp in [REF_GTF, EXON_STATUS_CSV, JUNCTION_CSV, GFF_FILE, FASTA_FILE]:
    if not os.path.exists(fp):
        raise FileNotFoundError(f"Missing file: {fp}")

exon_status = pd.read_csv(EXON_STATUS_CSV)
junction_status = pd.read_csv(JUNCTION_CSV)
fasta_lens = parse_fasta_lengths(FASTA_FILE)

if "in_ref" in exon_status.columns:
    exon_status["in_ref"] = exon_status["in_ref"].astype(str).str.lower().map(
        {"true": True, "false": False, "nan": False, "none": False}
    ).fillna(False)


#Load PacBio GFF exons
pb_rows = []
with open(GFF_FILE, "r") as fh:
    for line in fh:
        line = line.strip()
        if (not line) or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = parts
        if feature != "exon":
            continue
        ad = parse_gff_attributes(attrs)
        pb_rows.append({
            "seqnames": chrom,
            "source": source,
            "feature": feature,
            "start": int(start),
            "end": int(end),
            "strand": strand,
            "gene_id": ad.get("gene_id"),
            "transcript_id": ad.get("transcript_id"),
        })

pb_exons = pd.DataFrame(pb_rows)
if pb_exons.empty:
    raise ValueError("No exon rows found in selected_pbids.gff")
pb_exons["exon_key"] = pb_exons.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)


#Load reference ABCA4 exons
ref_rows = []
with gzip.open(REF_GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = parts
        if feature != "exon":
            continue
        ad = parse_gtf_attributes(attrs)
        if ad.get("gene_name") != "ABCA4":
            continue
        ref_rows.append({
            "seqnames": chrom,
            "start": int(start),
            "end": int(end),
            "strand": strand,
            "gene_id": ad.get("gene_id"),
            "gene_name": ad.get("gene_name"),
            "transcript_id": ad.get("transcript_id"),
            "transcript_name": ad.get("transcript_name"),
            "exon_number": ad.get("exon_number"),
        })

ref_exons = pd.DataFrame(ref_rows)
if ref_exons.empty:
    raise ValueError("No ABCA4 exons found in reference GTF")
ref_exons["exon_key"] = ref_exons.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)


#Precompute reference transcript structures
ref_tx_dict = {}
for tx, d in ref_exons.groupby("transcript_id"):
    strand = safe_mode(d["strand"].tolist())
    d_ord = ordered_exons(d, strand).copy()
    d_ord["tx_exon_index"] = range(1, d_ord.shape[0] + 1)

    #genomic-order version for intron calculations
    d_gen = d.sort_values(["start", "end"]).copy().reset_index(drop=True)
    introns = []
    for i in range(d_gen.shape[0] - 1):
        intron_start = int(d_gen.iloc[i]["end"]) + 1
        intron_end = int(d_gen.iloc[i + 1]["start"]) - 1
        if intron_start <= intron_end:
            introns.append({
                "intron_index_genomic": i + 1,
                "start": intron_start,
                "end": intron_end,
                "left_exon_end": int(d_gen.iloc[i]["end"]),
                "right_exon_start": int(d_gen.iloc[i + 1]["start"]),
            })

    ref_tx_dict[tx] = {
        "strand": strand,
        "transcript_name": safe_mode(d["transcript_name"].tolist()),
        "exons_tx_order": d_ord.reset_index(drop=True),
        "exons_genomic_order": d_gen,
        "exon_keys": set(d["exon_key"].tolist()),
        "introns": introns,
    }


#Merge exon/junction summaries
exon_status["exon_key"] = exon_status.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)
junction_status["junction_key"] = junction_status.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)

#optional cell type extraction
def get_cell_type_for_tx(tx):
    d = exon_status.loc[exon_status["transcript_id"] == tx].copy()
    cell_type = safe_mode(d["cell_type"].tolist()) if "cell_type" in d.columns else None
    if cell_type is None and ("track_id" in d.columns) and (not d.empty):
        vals = d["track_id"].dropna().astype(str).tolist()
        if vals:
            cell_type = vals[0].split("__")[0]
    return cell_type


#Choose best reference transcript for each PBID
def build_junctions_from_exons(df):
    d = df.sort_values(["start", "end"]).copy().reset_index(drop=True)
    rows = []
    for i in range(d.shape[0] - 1):
        jstart = int(d.iloc[i]["end"])
        jend = int(d.iloc[i + 1]["start"])
        rows.append(coord_key(jstart, jend))
    return rows


def choose_best_reference(pb_df):
    pb_exon_keys = set(pb_df["exon_key"].tolist())
    pb_junc_keys = set(build_junctions_from_exons(pb_df))

    best_tx = None
    best_score = None

    for ref_tx, info in ref_tx_dict.items():
        dref = info["exons_genomic_order"]
        ref_exon_keys = info["exon_keys"]
        ref_junc_keys = set(build_junctions_from_exons(dref))

        exact_exon_matches = len(pb_exon_keys & ref_exon_keys)
        exact_junction_matches = len(pb_junc_keys & ref_junc_keys)

        overlap_bp = 0
        for _, r1 in pb_df.iterrows():
            for _, r2 in dref.iterrows():
                overlap_bp += exon_overlap_len(r1["start"], r1["end"], r2["start"], r2["end"])

        score = (exact_junction_matches, exact_exon_matches, overlap_bp, -abs(pb_df.shape[0] - dref.shape[0]))
        if (best_score is None) or (score > best_score):
            best_score = score
            best_tx = ref_tx

    return best_tx, best_score


#Event extraction for SpliceAI comparison
def donor_acceptor_of_exon(exon_row, strand):
    """
    Return donor and acceptor genomic positions for this exon in transcript sense.
    '+' strand: acceptor=start, donor=end
    '-' strand: donor=start, acceptor=end
    """
    if strand == "+":
        return int(exon_row["end"]), int(exon_row["start"])   # donor, acceptor
    else:
        return int(exon_row["start"]), int(exon_row["end"])   # donor, acceptor


event_rows = []
summary_rows = []

for pbid, d_pb0 in pb_exons.groupby("transcript_id"):
    d_pb_gen = d_pb0.sort_values(["start", "end"]).copy().reset_index(drop=True)
    strand = safe_mode(d_pb_gen["strand"].tolist())
    d_pb_tx = ordered_exons(d_pb_gen, strand).copy()
    d_pb_tx["pb_exon_index"] = range(1, d_pb_tx.shape[0] + 1)

    cell_type = get_cell_type_for_tx(pbid)
    best_ref_tx, best_score = choose_best_reference(d_pb_gen)
    ref_info = ref_tx_dict[best_ref_tx]
    d_ref_tx = ref_info["exons_tx_order"].copy()
    d_ref_gen = ref_info["exons_genomic_order"].copy()

    fasta_len = fasta_lens.get(pbid)
    spliced_len = int((d_pb_gen["end"] - d_pb_gen["start"] + 1).sum())

    pb_exon_keys = set(d_pb_gen["exon_key"].tolist())
    ref_exon_keys = set(d_ref_gen["exon_key"].tolist())
    pb_junc_keys = set(build_junctions_from_exons(d_pb_gen))
    ref_junc_keys = set(build_junctions_from_exons(d_ref_gen))

    #read summary annotations
    d_ex_status = exon_status.loc[exon_status["transcript_id"] == pbid].copy()
    d_j_status = junction_status.loc[junction_status["transcript_id"] == pbid].copy()

    #Novel exons / pseudoexons
    novel_exon_count = 0
    if (not d_ex_status.empty) and ("diff_type" in d_ex_status.columns):
        novel_exon_count = int((d_ex_status["diff_type"] == "novel_exon").sum())

    for _, pex in d_pb_tx.iterrows():
        this_status = d_ex_status.loc[
            (d_ex_status["start"] == int(pex["start"])) &
            (d_ex_status["end"] == int(pex["end"]))
        ].copy()

        is_novel_exon = False
        if not this_status.empty and ("diff_type" in this_status.columns):
            is_novel_exon = (this_status["diff_type"] == "novel_exon").any()

        if is_novel_exon or (pex["exon_key"] not in ref_exon_keys):
            donor_pos, acceptor_pos = donor_acceptor_of_exon(pex, strand)
            event_rows.append({
                "transcript_id": pbid,
                "cell_type": cell_type,
                "reference_transcript_id": best_ref_tx,
                "reference_transcript_name": ref_info["transcript_name"],
                "strand": strand,
                "seqnames": pex["seqnames"],
                "event_type": "novel_exon",
                "event_subtype": "pseudoexon_or_unannotated_exon",
                "pb_exon_index": int(pex["pb_exon_index"]),
                "ref_exon_index": None,
                "observed_start": int(pex["start"]),
                "observed_end": int(pex["end"]),
                "observed_donor_pos": donor_pos,
                "observed_acceptor_pos": acceptor_pos,
                "reference_pos": None,
                "reference_donor_pos": None,
                "reference_acceptor_pos": None,
                "affected_region": f"{int(pex['start'])}-{int(pex['end'])}",
                "spliceai_expected_signal": "acceptor_gain_and_donor_gain",
                "note": "Novel exon observed in PBID; compare to predicted gained splice sites around exon boundaries."
            })

    #Alt donor/acceptor 
    for _, pex in d_pb_tx.iterrows():
        if pex["exon_key"] in ref_exon_keys:
            continue

        overlaps = d_ref_tx[
            (d_ref_tx["start"] <= int(pex["end"])) &
            (d_ref_tx["end"] >= int(pex["start"]))
        ].copy()

        if overlaps.empty:
            continue

        overlaps["ov"] = overlaps.apply(
            lambda r: exon_overlap_len(pex["start"], pex["end"], r["start"], r["end"]), axis=1
        )
        rbest = overlaps.sort_values("ov", ascending=False).iloc[0]

        pb_donor, pb_acceptor = donor_acceptor_of_exon(pex, strand)
        ref_donor, ref_acceptor = donor_acceptor_of_exon(rbest, strand)

        if pb_donor != ref_donor:
            event_rows.append({
                "transcript_id": pbid,
                "cell_type": cell_type,
                "reference_transcript_id": best_ref_tx,
                "reference_transcript_name": ref_info["transcript_name"],
                "strand": strand,
                "seqnames": pex["seqnames"],
                "event_type": "alt_donor",
                "event_subtype": "shifted_donor_site",
                "pb_exon_index": int(pex["pb_exon_index"]),
                "ref_exon_index": int(rbest["tx_exon_index"]),
                "observed_start": int(pex["start"]),
                "observed_end": int(pex["end"]),
                "observed_donor_pos": pb_donor,
                "observed_acceptor_pos": pb_acceptor,
                "reference_pos": ref_donor,
                "reference_donor_pos": ref_donor,
                "reference_acceptor_pos": ref_acceptor,
                "affected_region": f"exon_boundary_pb={int(pex['start'])}-{int(pex['end'])};ref={int(rbest['start'])}-{int(rbest['end'])}",
                "spliceai_expected_signal": "donor_gain_or_donor_loss",
                "note": "Observed donor differs from best overlapping reference exon."
            })

        if pb_acceptor != ref_acceptor:
            event_rows.append({
                "transcript_id": pbid,
                "cell_type": cell_type,
                "reference_transcript_id": best_ref_tx,
                "reference_transcript_name": ref_info["transcript_name"],
                "strand": strand,
                "seqnames": pex["seqnames"],
                "event_type": "alt_acceptor",
                "event_subtype": "shifted_acceptor_site",
                "pb_exon_index": int(pex["pb_exon_index"]),
                "ref_exon_index": int(rbest["tx_exon_index"]),
                "observed_start": int(pex["start"]),
                "observed_end": int(pex["end"]),
                "observed_donor_pos": pb_donor,
                "observed_acceptor_pos": pb_acceptor,
                "reference_pos": ref_acceptor,
                "reference_donor_pos": ref_donor,
                "reference_acceptor_pos": ref_acceptor,
                "affected_region": f"exon_boundary_pb={int(pex['start'])}-{int(pex['end'])};ref={int(rbest['start'])}-{int(rbest['end'])}",
                "spliceai_expected_signal": "acceptor_gain_or_acceptor_loss",
                "note": "Observed acceptor differs from best overlapping reference exon."
            })

    #Exon skipping
    d_ref_tx2 = d_ref_tx.copy().reset_index(drop=True)
    ref_left_by_end = {int(r["end"]): int(r["tx_exon_index"]) for _, r in d_ref_tx2.iterrows()}
    ref_right_by_start = {int(r["start"]): int(r["tx_exon_index"]) for _, r in d_ref_tx2.iterrows()}

    d_pb_gen2 = d_pb_gen.copy().reset_index(drop=True)
    for i in range(d_pb_gen2.shape[0] - 1):
        jstart = int(d_pb_gen2.iloc[i]["end"])
        jend = int(d_pb_gen2.iloc[i + 1]["start"])

        if (jstart in ref_left_by_end) and (jend in ref_right_by_start):
            left_idx = ref_left_by_end[jstart]
            right_idx = ref_right_by_start[jend]
            if abs(right_idx - left_idx) > 1:
                skipped = list(range(min(left_idx, right_idx) + 1, max(left_idx, right_idx)))
                event_rows.append({
                    "transcript_id": pbid,
                    "cell_type": cell_type,
                    "reference_transcript_id": best_ref_tx,
                    "reference_transcript_name": ref_info["transcript_name"],
                    "strand": strand,
                    "seqnames": safe_mode(d_pb_gen2["seqnames"].tolist()),
                    "event_type": "exon_skipping",
                    "event_subtype": "reference_internal_exon_skipping",
                    "pb_exon_index": None,
                    "ref_exon_index": ",".join(map(str, skipped)),
                    "observed_start": jstart,
                    "observed_end": jend,
                    "observed_donor_pos": jstart if strand == "+" else jend,
                    "observed_acceptor_pos": jend if strand == "+" else jstart,
                    "reference_pos": None,
                    "reference_donor_pos": None,
                    "reference_acceptor_pos": None,
                    "affected_region": f"skipped_ref_exons={','.join(map(str, skipped))}",
                    "spliceai_expected_signal": "donor_loss_or_acceptor_loss_or_cryptic_splicing",
                    "note": "PB junction joins non-adjacent reference exons."
                })

    #Intron retention / partial intron retention
    for _, pex in d_pb_tx.iterrows():
        for intr in ref_info["introns"]:
            intr_start = int(intr["start"])
            intr_end = int(intr["end"])
            if int(pex["start"]) < intr_start and int(pex["end"]) > intr_end:
                event_rows.append({
                    "transcript_id": pbid,
                    "cell_type": cell_type,
                    "reference_transcript_id": best_ref_tx,
                    "reference_transcript_name": ref_info["transcript_name"],
                    "strand": strand,
                    "seqnames": pex["seqnames"],
                    "event_type": "intron_retention",
                    "event_subtype": "full_bridge_across_reference_intron",
                    "pb_exon_index": int(pex["pb_exon_index"]),
                    "ref_exon_index": intr["intron_index_genomic"],
                    "observed_start": int(pex["start"]),
                    "observed_end": int(pex["end"]),
                    "observed_donor_pos": None,
                    "observed_acceptor_pos": None,
                    "reference_pos": f"{intr_start}-{intr_end}",
                    "reference_donor_pos": intr["left_exon_end"],
                    "reference_acceptor_pos": intr["right_exon_start"],
                    "affected_region": f"ref_intron={intr_start}-{intr_end}",
                    "spliceai_expected_signal": "donor_loss_or_acceptor_loss_or_intron_inclusion_related_signal",
                    "note": "PB exon spans across an entire reference intron."
                })
            elif intervals_overlap(pex["start"], pex["end"], intr_start, intr_end):
                event_rows.append({
                    "transcript_id": pbid,
                    "cell_type": cell_type,
                    "reference_transcript_id": best_ref_tx,
                    "reference_transcript_name": ref_info["transcript_name"],
                    "strand": strand,
                    "seqnames": pex["seqnames"],
                    "event_type": "partial_intron_retention",
                    "event_subtype": "partial_overlap_with_reference_intron",
                    "pb_exon_index": int(pex["pb_exon_index"]),
                    "ref_exon_index": intr["intron_index_genomic"],
                    "observed_start": int(pex["start"]),
                    "observed_end": int(pex["end"]),
                    "observed_donor_pos": None,
                    "observed_acceptor_pos": None,
                    "reference_pos": f"{intr_start}-{intr_end}",
                    "reference_donor_pos": intr["left_exon_end"],
                    "reference_acceptor_pos": intr["right_exon_start"],
                    "affected_region": f"ref_intron={intr_start}-{intr_end}",
                    "spliceai_expected_signal": "nearby_donor_or_acceptor_disruption_or_cryptic_site_gain",
                    "note": "PB exon partially overlaps a reference intron."
                })

    #Alternative first / last exon
    if (d_pb_tx.shape[0] > 0) and (d_ref_tx.shape[0] > 0):
        pb_first = d_pb_tx.iloc[0]
        pb_last = d_pb_tx.iloc[-1]
        ref_first = d_ref_tx.iloc[0]
        ref_last = d_ref_tx.iloc[-1]

        if (int(pb_first["start"]) != int(ref_first["start"])) or (int(pb_first["end"]) != int(ref_first["end"])):
            pb_donor, pb_acceptor = donor_acceptor_of_exon(pb_first, strand)
            ref_donor, ref_acceptor = donor_acceptor_of_exon(ref_first, strand)
            event_rows.append({
                "transcript_id": pbid,
                "cell_type": cell_type,
                "reference_transcript_id": best_ref_tx,
                "reference_transcript_name": ref_info["transcript_name"],
                "strand": strand,
                "seqnames": pb_first["seqnames"],
                "event_type": "alternative_first_exon",
                "event_subtype": "different_first_exon_or_boundary",
                "pb_exon_index": int(pb_first["pb_exon_index"]),
                "ref_exon_index": int(ref_first["tx_exon_index"]),
                "observed_start": int(pb_first["start"]),
                "observed_end": int(pb_first["end"]),
                "observed_donor_pos": pb_donor,
                "observed_acceptor_pos": pb_acceptor,
                "reference_pos": f"{int(ref_first['start'])}-{int(ref_first['end'])}",
                "reference_donor_pos": ref_donor,
                "reference_acceptor_pos": ref_acceptor,
                "affected_region": "first_exon",
                "spliceai_expected_signal": "possible_acceptor_or_donor_change_near_first_exon_boundary",
                "note": "PB first exon differs from best reference transcript."
            })

        if (int(pb_last["start"]) != int(ref_last["start"])) or (int(pb_last["end"]) != int(ref_last["end"])):
            pb_donor, pb_acceptor = donor_acceptor_of_exon(pb_last, strand)
            ref_donor, ref_acceptor = donor_acceptor_of_exon(ref_last, strand)
            event_rows.append({
                "transcript_id": pbid,
                "cell_type": cell_type,
                "reference_transcript_id": best_ref_tx,
                "reference_transcript_name": ref_info["transcript_name"],
                "strand": strand,
                "seqnames": pb_last["seqnames"],
                "event_type": "alternative_last_exon",
                "event_subtype": "different_last_exon_or_boundary",
                "pb_exon_index": int(pb_last["pb_exon_index"]),
                "ref_exon_index": int(ref_last["tx_exon_index"]),
                "observed_start": int(pb_last["start"]),
                "observed_end": int(pb_last["end"]),
                "observed_donor_pos": pb_donor,
                "observed_acceptor_pos": pb_acceptor,
                "reference_pos": f"{int(ref_last['start'])}-{int(ref_last['end'])}",
                "reference_donor_pos": ref_donor,
                "reference_acceptor_pos": ref_acceptor,
                "affected_region": "last_exon",
                "spliceai_expected_signal": "possible_acceptor_or_donor_change_near_last_exon_boundary",
                "note": "PB last exon differs from best reference transcript."
            })

    #transcript summary
    d_pb_gen2["tmp_exon_key"] = d_pb_gen2.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)
    novel_junction_count = 0
    if (not d_j_status.empty) and ("junction_category" in d_j_status.columns):
        novel_junction_count = int((d_j_status["junction_category"] == "novel").sum())

    event_types_for_tx = [r["event_type"] for r in event_rows if r["transcript_id"] == pbid]
    major_event = ";".join(sorted(set(event_types_for_tx))) if event_types_for_tx else "reference_like"

    summary_rows.append({
        "transcript_id": pbid,
        "cell_type": cell_type,
        "reference_transcript_id": best_ref_tx,
        "reference_transcript_name": ref_info["transcript_name"],
        "strand": strand,
        "seqnames": safe_mode(d_pb_gen["seqnames"].tolist()),
        "exon_count": d_pb_gen.shape[0],
        "junction_count": max(0, d_pb_gen.shape[0] - 1),
        "spliced_length_from_gff": spliced_len,
        "fasta_length": fasta_len,
        "length_match_fasta_vs_spliced": None if fasta_len is None else (int(fasta_len) == int(spliced_len)),
        "exact_exon_match_count_to_ref": len(pb_exon_keys & ref_exon_keys),
        "exact_junction_match_count_to_ref": len(pb_junc_keys & ref_junc_keys),
        "novel_exon_count_from_summary": novel_exon_count,
        "novel_junction_count_from_summary": novel_junction_count,
        "event_count_for_spliceai": sum(1 for r in event_rows if r["transcript_id"] == pbid),
        "major_event_for_spliceai": major_event,
    })


#Save outputs
events_df = pd.DataFrame(event_rows)
summary_df = pd.DataFrame(summary_rows)

if not events_df.empty:
    events_df = events_df.sort_values(
        ["transcript_id", "event_type", "pb_exon_index", "observed_start", "observed_end"],
        na_position="last"
    ).reset_index(drop=True)

summary_df = summary_df.sort_values(
    ["major_event_for_spliceai", "cell_type", "transcript_id"]
).reset_index(drop=True)

events_file = os.path.join(OUT_DIR, "abca4_spliceai_comparison_events.tsv")
summary_file = os.path.join(OUT_DIR, "abca4_spliceai_transcript_summary.tsv")
report_file = os.path.join(OUT_DIR, "abca4_spliceai_ready_report.txt")

events_df.to_csv(events_file, sep="\t", index=False)
summary_df.to_csv(summary_file, sep="\t", index=False)

with open(report_file, "w") as out:
    out.write("ABCA4 SpliceAI-ready event table\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"REF_GTF = {REF_GTF}\n")
    out.write(f"EXON_STATUS_CSV = {EXON_STATUS_CSV}\n")
    out.write(f"JUNCTION_CSV = {JUNCTION_CSV}\n")
    out.write(f"GFF_FILE = {GFF_FILE}\n")
    out.write(f"FASTA_FILE = {FASTA_FILE}\n\n")
    out.write(f"Transcript summary rows: {summary_df.shape[0]}\n")
    out.write(f"Event rows: {events_df.shape[0]}\n\n")
    if not events_df.empty:
        out.write("Event type counts:\n")
        for k, v in events_df["event_type"].value_counts().items():
            out.write(f"{k}\t{v}\n")

print("Done.")
print("Saved:")
print(events_file)
print(summary_file)
print(report_file)

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################