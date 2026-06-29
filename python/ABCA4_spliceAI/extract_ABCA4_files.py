########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import defaultdict

#paths
work_dir = Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI")
fasta_file = Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/orfanage/RetinaMerge.filtered_lite.transcripts.fa")
gff_file = Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/original_from_Ray/ray_pigeon_output/classify.retina.merge.GRCh38.p14.with.cage.polya.FL/RetinaMerge.sorted.filtered_lite.gff")
target_prefix = "PB.5877."
out_dir = work_dir / "extracted_pbids"
per_fasta_dir = out_dir / "per_pbid_fasta"
per_gff_dir = out_dir / "per_pbid_gff"
out_dir.mkdir(parents=True, exist_ok=True)
per_fasta_dir.mkdir(parents=True, exist_ok=True)
per_gff_dir.mkdir(parents=True, exist_ok=True)

#define functions
def parse_fasta(path):
    """
    Minimal FASTA parser.
    Uses the first token after '>' as record id.
    """
    seqs = {}
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            seqs[header] = "".join(seq_chunks)
    return seqs

def wrap_seq(seq, width=60):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def line_matches_pbid(line, pbid):
    return pbid in line

#Read FASTA and determine all matching PBids
print("Reading FASTA ...")
all_seqs = parse_fasta(fasta_file)

target_pbids = sorted([pbid for pbid in all_seqs if pbid.startswith(target_prefix)])
pbid_set = set(target_pbids)

print(f"Found {len(target_pbids)} PBids with prefix {target_prefix}")

if len(target_pbids) == 0:
    raise ValueError(f"No PBids found with prefix: {target_prefix}")

#Write PBid list
pbid_txt = out_dir / "pbids.txt"
with open(pbid_txt, "w") as f:
    for pbid in target_pbids:
        f.write(pbid + "\n")

print(f"Wrote PBid list: {pbid_txt}")

#Extract FASTA
selected_fasta = out_dir / "selected_pbids.fa"
found_fasta = []
missing_fasta = []

with open(selected_fasta, "w") as fout:
    for pbid in target_pbids:
        if pbid in all_seqs:
            seq = all_seqs[pbid]
            fout.write(f">{pbid}\n")
            fout.write(wrap_seq(seq) + "\n")
            found_fasta.append(pbid)

            per_file = per_fasta_dir / f"{pbid}.fa"
            with open(per_file, "w") as pf:
                pf.write(f">{pbid}\n")
                pf.write(wrap_seq(seq) + "\n")
        else:
            missing_fasta.append(pbid)

print(f"FASTA found: {len(found_fasta)} / {len(target_pbids)}")
if missing_fasta:
    print("FASTA missing PBids:")
    for x in missing_fasta:
        print("   ", x)

#Extract GFF
print("Reading GFF and extracting matching lines ...")

selected_gff = out_dir / "selected_pbids.gff"
gff_lines_by_pbid = defaultdict(list)
header_lines = []

with open(gff_file) as fin, open(selected_gff, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            header_lines.append(line)
            fout.write(line)
            continue

        matched_any = False
        for pbid in target_pbids:
            if line_matches_pbid(line, pbid):
                matched_any = True
                gff_lines_by_pbid[pbid].append(line)

        if matched_any:
            fout.write(line)

#write per-PBid gff
missing_gff = []
for pbid in target_pbids:
    per_file = per_gff_dir / f"{pbid}.gff"
    lines = gff_lines_by_pbid.get(pbid, [])
    if len(lines) == 0:
        missing_gff.append(pbid)

    with open(per_file, "w") as f:
        for h in header_lines:
            f.write(h)
        for x in lines:
            f.write(x)

print(f"GFF found: {len(target_pbids) - len(missing_gff)} / {len(target_pbids)}")
if missing_gff:
    print("GFF missing PBids:")
    for x in missing_gff:
        print("   ", x)

#Summary table
summary_file = out_dir / "summary.tsv"
with open(summary_file, "w") as f:
    f.write("PBid\tfound_in_fasta\tfound_in_gff\tfasta_length\tgff_line_count\n")
    for pbid in target_pbids:
        fasta_ok = "yes" if pbid in all_seqs else "no"
        gff_ok = "yes" if len(gff_lines_by_pbid.get(pbid, [])) > 0 else "no"
        fasta_len = len(all_seqs[pbid]) if pbid in all_seqs else 0
        gff_n = len(gff_lines_by_pbid.get(pbid, []))
        f.write(f"{pbid}\t{fasta_ok}\t{gff_ok}\t{fasta_len}\t{gff_n}\n")

print(f"Wrote summary: {summary_file}")
print("Done.")

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################