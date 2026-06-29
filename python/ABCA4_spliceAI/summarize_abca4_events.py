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
FASTA_FILE = os.path.join(EXTRACTED_DIR, "selected_pbids.fa")
GFF_FILE = os.path.join(EXTRACTED_DIR, "selected_pbids.gff")

OUT_DIR = os.path.join("/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI", "ABCA4_event_summary_v44")
os.makedirs(OUT_DIR, exist_ok=True)

#define functions
def parse_gtf_attributes(attr_str):
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


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


def parse_gff_attributes(attr_str):
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


def coord_key(start, end):
    return f"{int(start)}-{int(end)}"


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


def build_exon_chain(df, strand):
    d = ordered_exons(df, strand)
    return "|".join([coord_key(s, e) for s, e in zip(d["start"], d["end"])])


def build_junction_df_from_exons(df, strand, transcript_id):
    d = ordered_exons(df, strand)
    rows = []
    if d.shape[0] < 2:
        return pd.DataFrame(columns=["transcript_id", "start", "end", "junction_key", "strand"])
    for i in range(d.shape[0] - 1):
        e1s, e1e = int(d.loc[i, "start"]), int(d.loc[i, "end"])
        e2s, e2e = int(d.loc[i + 1, "start"]), int(d.loc[i + 1, "end"])
        # always keep genomic coords as lower exon end -> next exon start in genomic order
        # so use genomic order for coordinates
        g = df.sort_values(["start", "end"]).copy()
    g = df.sort_values(["start", "end"]).copy()
    rows = []
    for i in range(g.shape[0] - 1):
        jstart = int(g.iloc[i]["end"])
        jend = int(g.iloc[i + 1]["start"])
        rows.append({
            "transcript_id": transcript_id,
            "start": jstart,
            "end": jend,
            "junction_key": coord_key(jstart, jend),
            "strand": strand
        })
    return pd.DataFrame(rows)


def exon_overlap_len(a_start, a_end, b_start, b_end):
    x1 = max(int(a_start), int(b_start))
    x2 = min(int(a_end), int(b_end))
    return max(0, x2 - x1 + 1)


def intervals_overlap(a_start, a_end, b_start, b_end):
    return exon_overlap_len(a_start, a_end, b_start, b_end) > 0


def load_pacbio_gff_exons(gff_file):
    exons = []
    transcripts = []
    with open(gff_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if (not line) or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            start = int(start)
            end = int(end)
            ad = parse_gff_attributes(attrs)
            row = {
                "seqnames": chrom,
                "source": source,
                "feature": feature,
                "start": start,
                "end": end,
                "score": score,
                "strand": strand,
                "frame": frame,
                "gene_id": ad.get("gene_id"),
                "transcript_id": ad.get("transcript_id"),
            }
            if feature == "exon":
                exons.append(row)
            elif feature == "transcript":
                transcripts.append(row)
    exon_df = pd.DataFrame(exons)
    tx_df = pd.DataFrame(transcripts)
    if exon_df.empty:
        raise ValueError("No exon entries found in selected_pbids.gff")
    exon_df["exon_key"] = exon_df.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)
    return exon_df, tx_df


def load_reference_abca4_from_gtf(gtf_gz):
    rows = []
    open_func = gzip.open if gtf_gz.endswith(".gz") else open
    with open_func(gtf_gz, "rt") as fh:
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
            gene_name = ad.get("gene_name")
            if gene_name != "ABCA4":
                continue
            rows.append({
                "seqnames": chrom,
                "source": source,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "score": score,
                "strand": strand,
                "frame": frame,
                "gene_id": ad.get("gene_id"),
                "gene_name": gene_name,
                "transcript_id": ad.get("transcript_id"),
                "transcript_name": ad.get("transcript_name"),
                "exon_number": ad.get("exon_number"),
                "exon_id": ad.get("exon_id"),
            })
    ref_exons = pd.DataFrame(rows)
    if ref_exons.empty:
        raise ValueError("No ABCA4 exon entries found in reference GTF.")
    ref_exons["exon_key"] = ref_exons.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)
    return ref_exons


def assign_global_reference_exon_index(ref_exons):
    """
    Build a nonredundant ABCA4 reference exon catalog across all Gencode transcripts.
    """
    d = ref_exons[["seqnames", "start", "end", "strand", "exon_key"]].drop_duplicates().copy()
    strand = safe_mode(d["strand"].tolist())
    d = d.sort_values(["start", "end"]).copy()
    d["global_ref_exon_index_genomic"] = range(1, d.shape[0] + 1)
    if strand == "-":
        d["global_ref_exon_index_tx"] = list(range(d.shape[0], 0, -1))
    else:
        d["global_ref_exon_index_tx"] = range(1, d.shape[0] + 1)
    return d


def transcript_junction_chain_from_exons(df, strand):
    d = ordered_exons(df, strand)
    if d.shape[0] < 2:
        return ""
    coords = []
    g = df.sort_values(["start", "end"]).copy()
    for i in range(g.shape[0] - 1):
        coords.append(coord_key(int(g.iloc[i]["end"]), int(g.iloc[i + 1]["start"])))
    if strand == "-":
        coords = coords[::-1]
    return "|".join(coords)


def choose_best_reference_for_pbid(pb_exons, ref_exons):
    """
    Pick the best ABCA4 reference transcript by:
    1) number of exact exon matches
    2) number of exact junction matches
    3) total overlapped exonic bp
    """
    strand = safe_mode(pb_exons["strand"].tolist())
    pb_tx = safe_mode(pb_exons["transcript_id"].tolist())
    pb_exon_keys = set(pb_exons["exon_key"].tolist())

    pb_junc_df = build_junction_df_from_exons(pb_exons, strand, pb_tx)
    pb_junc_keys = set(pb_junc_df["junction_key"].tolist()) if not pb_junc_df.empty else set()

    best = None
    best_score = None

    for ref_tx, dref in ref_exons.groupby("transcript_id"):
        ref_exon_keys = set(dref["exon_key"].tolist())

        ref_junc_df = build_junction_df_from_exons(dref, strand, ref_tx)
        ref_junc_keys = set(ref_junc_df["junction_key"].tolist()) if not ref_junc_df.empty else set()

        exact_exon_matches = len(pb_exon_keys & ref_exon_keys)
        exact_junction_matches = len(pb_junc_keys & ref_junc_keys)

        overlap_bp = 0
        for _, r1 in pb_exons.iterrows():
            for _, r2 in dref.iterrows():
                overlap_bp += exon_overlap_len(r1["start"], r1["end"], r2["start"], r2["end"])

        score = (exact_junction_matches, exact_exon_matches, overlap_bp, -abs(pb_exons.shape[0] - dref.shape[0]))

        if (best_score is None) or (score > best_score):
            best_score = score
            best = {
                "reference_transcript_id": ref_tx,
                "reference_transcript_name": safe_mode(dref["transcript_name"].tolist()),
                "exact_junction_matches": exact_junction_matches,
                "exact_exon_matches": exact_exon_matches,
                "overlap_bp": overlap_bp,
                "ref_exon_count": dref.shape[0],
                "score_tuple": score,
            }

    return best


def classify_against_reference(pb_exons, ref_exons_one_tx, global_ref_catalog, exon_status_tx=None, junction_status_tx=None):
    """
    Event calling against best matched reference transcript.
    """
    strand = safe_mode(pb_exons["strand"].tolist())
    pb_ord = ordered_exons(pb_exons, strand)
    ref_ord = ordered_exons(ref_exons_one_tx, strand)

    pb_exon_keys = set(pb_ord["exon_key"].tolist())
    ref_exon_keys = set(ref_ord["exon_key"].tolist())

    pb_junc_df = build_junction_df_from_exons(pb_exons, strand, safe_mode(pb_exons["transcript_id"].tolist()))
    ref_junc_df = build_junction_df_from_exons(ref_exons_one_tx, strand, safe_mode(ref_exons_one_tx["transcript_id"].tolist()))

    pb_junc_keys = set(pb_junc_df["junction_key"].tolist()) if not pb_junc_df.empty else set()
    ref_junc_keys = set(ref_junc_df["junction_key"].tolist()) if not ref_junc_df.empty else set()

    exact_reference_match = (pb_exon_keys == ref_exon_keys)

    #basic novel counts from provided summary files
    novel_exon_count = 0
    known_exon_count = 0
    if exon_status_tx is not None and (not exon_status_tx.empty):
        novel_exon_count = int((exon_status_tx["diff_type"] == "novel_exon").sum())
        known_exon_count = int((exon_status_tx["diff_type"] == "known_exon").sum())

    novel_junction_count = 0
    known_junction_count = 0
    if junction_status_tx is not None and (not junction_status_tx.empty):
        novel_junction_count = int((junction_status_tx["junction_category"] == "novel").sum())
        known_junction_count = int((junction_status_tx["junction_category"] == "known").sum())

    #terminal differences
    alt_first_exon = False
    alt_last_exon = False
    if pb_ord.shape[0] > 0 and ref_ord.shape[0] > 0:
        pb_first = pb_ord.iloc[0]
        pb_last = pb_ord.iloc[-1]
        ref_first = ref_ord.iloc[0]
        ref_last = ref_ord.iloc[-1]
        if (int(pb_first["start"]) != int(ref_first["start"])) or (int(pb_first["end"]) != int(ref_first["end"])):
            alt_first_exon = True
        if (int(pb_last["start"]) != int(ref_last["start"])) or (int(pb_last["end"]) != int(ref_last["end"])):
            alt_last_exon = True

    #exon skipping:
    #if a PB junction connects two reference exons that are not adjacent in the matched reference transcript
    exon_skip_count = 0
    ref_exon_order = ordered_exons(ref_exons_one_tx, strand).reset_index(drop=True).copy()
    ref_exon_order["ref_order"] = range(1, ref_exon_order.shape[0] + 1)

    for _, pj in pb_junc_df.iterrows():
        left_ref = ref_exon_order.loc[ref_exon_order["end"] == int(pj["start"])]
        right_ref = ref_exon_order.loc[ref_exon_order["start"] == int(pj["end"])]
        if left_ref.shape[0] == 1 and right_ref.shape[0] == 1:
            lidx = int(left_ref["ref_order"].iloc[0])
            ridx = int(right_ref["ref_order"].iloc[0])
            if abs(ridx - lidx) > 1:
                exon_skip_count += abs(ridx - lidx) - 1

    has_exon_skipping = exon_skip_count > 0

    #alt donor / acceptor:
    #compare unmatched PB exon boundaries to overlapping reference exons
    alt_donor = False
    alt_acceptor = False

    for _, pex in pb_ord.iterrows():
        if pex["exon_key"] in ref_exon_keys:
            continue
        overlaps = ref_ord[
            (ref_ord["start"] <= int(pex["end"])) &
            (ref_ord["end"] >= int(pex["start"]))
        ].copy()
        if overlaps.empty:
            continue
        #use highest overlap reference exon
        overlaps["ov"] = overlaps.apply(lambda r: exon_overlap_len(pex["start"], pex["end"], r["start"], r["end"]), axis=1)
        rbest = overlaps.sort_values("ov", ascending=False).iloc[0]

        if strand == "+":
            #donor = exon end, acceptor = exon start
            if int(pex["end"]) != int(rbest["end"]):
                alt_donor = True
            if int(pex["start"]) != int(rbest["start"]):
                alt_acceptor = True
        else:
            #on negative strand donor/acceptor are swapped in genomic boundary sense
            #transcript-aware:
            #donor = exon start, acceptor = exon end
            if int(pex["start"]) != int(rbest["start"]):
                alt_donor = True
            if int(pex["end"]) != int(rbest["end"]):
                alt_acceptor = True

    #intron retention / partial IR tendency:
    #if one PB exon spans across boundaries of >=2 reference exons or bridges a known reference intron
    intron_retention = False
    partial_intron_retention = False

    ref_ord2 = ref_ord.sort_values(["start", "end"]).copy().reset_index(drop=True)
    for _, pex in pb_ord.iterrows():
        overlapping_ref_exons = ref_ord2[
            (ref_ord2["start"] <= int(pex["end"])) &
            (ref_ord2["end"] >= int(pex["start"]))
        ].copy()

        if overlapping_ref_exons.shape[0] >= 2:
            intron_retention = True

        #bridge intron between adjacent reference exons
        for i in range(ref_ord2.shape[0] - 1):
            intron_start = int(ref_ord2.iloc[i]["end"]) + 1
            intron_end = int(ref_ord2.iloc[i + 1]["start"]) - 1
            if intron_start <= intron_end:
                if int(pex["start"]) < intron_start and int(pex["end"]) > intron_end:
                    intron_retention = True
                elif intervals_overlap(pex["start"], pex["end"], intron_start, intron_end):
                    partial_intron_retention = True

    #novel exon presence from exon_status
    has_novel_exon = novel_exon_count > 0
    has_novel_junction = novel_junction_count > 0

    #assign major event conservatively
    if exact_reference_match and (not has_novel_exon) and (not has_novel_junction):
        major_event = "reference_match"
    else:
        event_tags = []
        if has_novel_exon:
            event_tags.append("novel_exon")
        if has_novel_junction:
            event_tags.append("novel_junction")
        if has_exon_skipping:
            event_tags.append("exon_skipping")
        if alt_donor and alt_acceptor:
            event_tags.append("alt_donor_acceptor")
        elif alt_donor:
            event_tags.append("alt_donor")
        elif alt_acceptor:
            event_tags.append("alt_acceptor")
        if intron_retention:
            event_tags.append("intron_retention")
        elif partial_intron_retention:
            event_tags.append("partial_intron_retention")
        if alt_first_exon:
            event_tags.append("alternative_first_exon")
        if alt_last_exon:
            event_tags.append("alternative_last_exon")

        if not event_tags:
            major_event = "reference_like_but_shifted"
        else:
            major_event = ";".join(event_tags)

    return {
        "exact_reference_match": exact_reference_match,
        "novel_exon_count": novel_exon_count,
        "known_exon_count": known_exon_count,
        "novel_junction_count": novel_junction_count,
        "known_junction_count": known_junction_count,
        "has_novel_exon": has_novel_exon,
        "has_novel_junction": has_novel_junction,
        "has_exon_skipping": has_exon_skipping,
        "exon_skip_count_estimate": exon_skip_count,
        "alt_donor": alt_donor,
        "alt_acceptor": alt_acceptor,
        "intron_retention": intron_retention,
        "partial_intron_retention": partial_intron_retention,
        "alternative_first_exon": alt_first_exon,
        "alternative_last_exon": alt_last_exon,
        "major_event": major_event,
        "pb_exon_chain": build_exon_chain(pb_exons, strand),
        "ref_exon_chain": build_exon_chain(ref_exons_one_tx, strand),
        "pb_junction_chain": transcript_junction_chain_from_exons(pb_exons, strand),
        "ref_junction_chain": transcript_junction_chain_from_exons(ref_exons_one_tx, strand),
    }


#Load files
for fp in [REF_GTF, EXON_STATUS_CSV, JUNCTION_CSV, FASTA_FILE, GFF_FILE]:
    if not os.path.exists(fp):
        raise FileNotFoundError(f"Missing file: {fp}")

print("Loading inputs ...")

exon_status = pd.read_csv(EXON_STATUS_CSV)
junction_status = pd.read_csv(JUNCTION_CSV)
pb_exons, pb_txs = load_pacbio_gff_exons(GFF_FILE)
fasta_lens = parse_fasta_lengths(FASTA_FILE)
ref_abca4 = load_reference_abca4_from_gtf(REF_GTF)
global_ref_catalog = assign_global_reference_exon_index(ref_abca4)

# normalize boolean-ish in_ref
if "in_ref" in exon_status.columns:
    exon_status["in_ref"] = exon_status["in_ref"].astype(str).str.lower().map({
        "true": True, "false": False, "nan": False, "none": False
    }).fillna(False)

exon_status["exon_key"] = exon_status.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)
junction_status["junction_key"] = junction_status.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)
ref_abca4["exon_key"] = ref_abca4.apply(lambda r: coord_key(r["start"], r["end"]), axis=1)


#Annotate PB exons with global reference exon index if exact coord match
global_ref_map = global_ref_catalog.set_index("exon_key")["global_ref_exon_index_tx"].to_dict()
pb_exons["global_ref_exon_index_tx"] = pb_exons["exon_key"].map(global_ref_map)

exon_status_sub = exon_status[[
    "transcript_id", "seqnames", "start", "end", "strand",
    "cell_type", "in_ref", "diff_type", "track_id", "transcript_to_plot"
]].drop_duplicates()

pb_exons_annot = pb_exons.merge(
    exon_status_sub,
    how="left",
    on=["transcript_id", "seqnames", "start", "end", "strand"]
)

pb_exons_annot["diff_type"] = pb_exons_annot["diff_type"].fillna("unannotated_in_summary")
pb_exons_annot["in_ref"] = pb_exons_annot["in_ref"].fillna(False)


#Per-PBID summary
summary_rows = []
detail_rows = []

pbids = sorted(pb_exons["transcript_id"].dropna().unique().tolist())

for tx in pbids:
    d_pb = pb_exons_annot.loc[pb_exons_annot["transcript_id"] == tx].copy()
    if d_pb.empty:
        continue

    strand = safe_mode(d_pb["strand"].tolist())
    cell_type = safe_mode(d_pb["cell_type"].tolist())
    if cell_type is None:
        raw = exon_status.loc[exon_status["transcript_id"] == tx]
        if not raw.empty and "track_id" in raw.columns:
            vals = raw["track_id"].dropna().astype(str).tolist()
            if vals:
                cell_type = vals[0].split("__")[0]

    best_ref = choose_best_reference_for_pbid(d_pb, ref_abca4)
    best_ref_tx = best_ref["reference_transcript_id"]
    d_ref = ref_abca4.loc[ref_abca4["transcript_id"] == best_ref_tx].copy()

    d_ex_status_tx = exon_status.loc[exon_status["transcript_id"] == tx].copy()
    d_junc_status_tx = junction_status.loc[junction_status["transcript_id"] == tx].copy()

    cls = classify_against_reference(
        pb_exons=d_pb,
        ref_exons_one_tx=d_ref,
        global_ref_catalog=global_ref_catalog,
        exon_status_tx=d_ex_status_tx,
        junction_status_tx=d_junc_status_tx
    )

    exon_count = d_pb.shape[0]
    junction_count = max(0, exon_count - 1)
    tx_start = int(d_pb["start"].min())
    tx_end = int(d_pb["end"].max())
    spliced_len = int((d_pb["end"] - d_pb["start"] + 1).sum())
    fasta_len = fasta_lens.get(tx)
    length_match = None if fasta_len is None else (int(fasta_len) == int(spliced_len))

    priority = "medium"
    if cls["has_novel_exon"] and cls["has_novel_junction"]:
        priority = "high"
    elif cls["major_event"] == "reference_match":
        priority = "low"
    elif cls["intron_retention"] or cls["has_exon_skipping"]:
        priority = "high"
    elif cls["alt_donor"] or cls["alt_acceptor"]:
        priority = "high"

    summary_rows.append({
        "transcript_id": tx,
        "cell_type": cell_type,
        "strand": strand,
        "seqnames": safe_mode(d_pb["seqnames"].tolist()),
        "tx_start": tx_start,
        "tx_end": tx_end,
        "genomic_span_bp": tx_end - tx_start + 1,
        "spliced_length_from_gff": spliced_len,
        "fasta_length": fasta_len,
        "length_match_fasta_vs_spliced": length_match,
        "exon_count": exon_count,
        "junction_count": junction_count,

        "best_reference_transcript_id": best_ref["reference_transcript_id"],
        "best_reference_transcript_name": best_ref["reference_transcript_name"],
        "best_reference_exon_count": best_ref["ref_exon_count"],
        "exact_exon_matches_to_best_ref": best_ref["exact_exon_matches"],
        "exact_junction_matches_to_best_ref": best_ref["exact_junction_matches"],
        "overlap_bp_to_best_ref": best_ref["overlap_bp"],

        "exact_reference_match": cls["exact_reference_match"],
        "novel_exon_count": cls["novel_exon_count"],
        "known_exon_count": cls["known_exon_count"],
        "novel_junction_count": cls["novel_junction_count"],
        "known_junction_count": cls["known_junction_count"],
        "has_novel_exon": cls["has_novel_exon"],
        "has_novel_junction": cls["has_novel_junction"],
        "has_exon_skipping": cls["has_exon_skipping"],
        "exon_skip_count_estimate": cls["exon_skip_count_estimate"],
        "alt_donor": cls["alt_donor"],
        "alt_acceptor": cls["alt_acceptor"],
        "intron_retention": cls["intron_retention"],
        "partial_intron_retention": cls["partial_intron_retention"],
        "alternative_first_exon": cls["alternative_first_exon"],
        "alternative_last_exon": cls["alternative_last_exon"],
        "major_event": cls["major_event"],
        "priority_for_followup": priority,

        "pb_exon_chain": cls["pb_exon_chain"],
        "ref_exon_chain": cls["ref_exon_chain"],
        "pb_junction_chain": cls["pb_junction_chain"],
        "ref_junction_chain": cls["ref_junction_chain"],
    })

    detail_rows.append({
        "transcript_id": tx,
        "best_reference_transcript_id": best_ref["reference_transcript_id"],
        "best_reference_transcript_name": best_ref["reference_transcript_name"],
        "score_tuple": str(best_ref["score_tuple"]),
        "pb_exon_chain": cls["pb_exon_chain"],
        "ref_exon_chain": cls["ref_exon_chain"],
        "pb_junction_chain": cls["pb_junction_chain"],
        "ref_junction_chain": cls["ref_junction_chain"],
        "major_event": cls["major_event"],
    })

summary_df = pd.DataFrame(summary_rows)
detail_df = pd.DataFrame(detail_rows)

summary_df = summary_df.sort_values(
    ["priority_for_followup", "major_event", "cell_type", "transcript_id"],
    ascending=[True, True, True, True]
).reset_index(drop=True)


#Save outputs
summary_file = os.path.join(OUT_DIR, "abca4_isoform_event_summary.tsv")
detail_file = os.path.join(OUT_DIR, "abca4_transcript_vs_reference_details.tsv")
ref_exon_file = os.path.join(OUT_DIR, "abca4_reference_exons_v44.tsv")
pb_exon_annot_file = os.path.join(OUT_DIR, "abca4_exon_table_annotated.tsv")
report_file = os.path.join(OUT_DIR, "abca4_run_report.txt")

summary_df.to_csv(summary_file, sep="\t", index=False)
detail_df.to_csv(detail_file, sep="\t", index=False)
global_ref_catalog.to_csv(ref_exon_file, sep="\t", index=False)
pb_exons_annot.to_csv(pb_exon_annot_file, sep="\t", index=False)

with open(report_file, "w") as out:
    out.write("ABCA4 isoform event summary report\n")
    out.write("=" * 60 + "\n\n")
    out.write("Input files:\n")
    out.write(f"REF_GTF: {REF_GTF}\n")
    out.write(f"EXON_STATUS_CSV: {EXON_STATUS_CSV}\n")
    out.write(f"JUNCTION_CSV: {JUNCTION_CSV}\n")
    out.write(f"FASTA_FILE: {FASTA_FILE}\n")
    out.write(f"GFF_FILE: {GFF_FILE}\n\n")

    out.write("Output files:\n")
    out.write(f"{summary_file}\n")
    out.write(f"{detail_file}\n")
    out.write(f"{ref_exon_file}\n")
    out.write(f"{pb_exon_annot_file}\n\n")

    out.write("Counts:\n")
    out.write(f"PBIDs analyzed: {summary_df.shape[0]}\n")
    out.write(f"Reference ABCA4 transcripts in v44: {ref_abca4['transcript_id'].nunique()}\n")
    out.write(f"Reference ABCA4 unique exons: {global_ref_catalog.shape[0]}\n\n")

    out.write("Major event counts:\n")
    for k, v in summary_df["major_event"].value_counts().items():
        out.write(f"{k}\t{v}\n")
    out.write("\n")

print("Done.")
print("Saved:")
print(summary_file)
print(detail_file)
print(ref_exon_file)
print(pb_exon_annot_file)
print(report_file)

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################