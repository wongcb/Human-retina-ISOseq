########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pandas as pd

#paths
INPUT_CSV = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/clinVar_ABCA4.csv"

REFERENCE_FASTA = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/GRCh38p14/GRCh38.p14.genome.fa"

OUT_DIR = "/mnt/d/bioinformatic/isoseq/processed/from_Ray/revision/splicAI/clinvar_spliceai_input"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_VCF = os.path.join(OUT_DIR, "clinvar_ABCA4.spliceai_input.vcf")
OUT_EXCLUDED = os.path.join(OUT_DIR, "clinvar_ABCA4.excluded.tsv")
OUT_INCLUDED = os.path.join(OUT_DIR, "clinvar_ABCA4.included.tsv")
OUT_REPORT = os.path.join(OUT_DIR, "clinvar_ABCA4.vcf_build_report.txt")


#define functions
class FastaReader:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.seqs = {}
        self._load()

    def _load(self):
        chrom = None
        seq_parts = []
        with open(self.fasta_path, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if chrom is not None:
                        self.seqs[chrom] = "".join(seq_parts).upper()
                    chrom = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())
            if chrom is not None:
                self.seqs[chrom] = "".join(seq_parts).upper()

    def fetch_base_1based(self, chrom, pos1):
        if chrom not in self.seqs:
            raise KeyError(f"Chromosome not found in FASTA: {chrom}")
        seq = self.seqs[chrom]
        if pos1 < 1 or pos1 > len(seq):
            raise ValueError(f"Position out of range: {chrom}:{pos1}")
        return seq[pos1 - 1]

    def fetch_seq_1based_closed(self, chrom, start1, end1):
        if chrom not in self.seqs:
            raise KeyError(f"Chromosome not found in FASTA: {chrom}")
        seq = self.seqs[chrom]
        if start1 < 1 or end1 > len(seq) or start1 > end1:
            raise ValueError(f"Interval out of range: {chrom}:{start1}-{end1}")
        return seq[start1 - 1:end1]


#Parsers
def parse_spdi(spdi):
    """
    SPDI format:
    Sequence:Position:Deleted:Inserted

    Position is 0-based interbase.
    """
    if pd.isna(spdi):
        return None
    s = str(spdi).strip()
    m = re.match(r"^([^:]+):(\d+):([^:]*):([^:]*)$", s)
    if not m:
        return None
    return {
        "seq_id": m.group(1),
        "pos0": int(m.group(2)),
        "deleted": m.group(3),
        "inserted": m.group(4),
    }


def chrom_to_fasta_name(chrom_value):
    """
    ClinVar GRCh38Chromosome uses 1,2,...,X,Y,MT
    Our FASTA uses chr1, chr2, ...
    """
    if pd.isna(chrom_value):
        return None
    c = str(chrom_value).strip()
    if c.startswith("chr"):
        return c
    if c == "MT":
        return "chrM"
    return "chr" + c


def trim_common(ref, alt):
    """
    Minimal trimming of shared prefix/suffix.
    Return trimmed ref, alt, left_shift
    """
    left_shift = 0

    #trim common prefix
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        left_shift += 1

    #trim common suffix
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    return ref, alt, left_shift


def spdi_to_vcf_record(spdi_obj, chrom_fasta, fasta):
    """
    Convert SPDI to left-anchored VCF alleles.
    Returns dict with CHROM, POS, REF, ALT
    or raises exception if conversion fails.
    """
    pos0 = spdi_obj["pos0"]
    deleted = spdi_obj["deleted"].upper()
    inserted = spdi_obj["inserted"].upper()

    #SNV / MNV / delins before normalization
    ref_trim, alt_trim, left_shift = trim_common(deleted, inserted)

    #SPDI deleted sequence starts at 1-based pos0+1
    deleted_start1 = pos0 + 1 + left_shift

    #substitution after trimming
    if len(ref_trim) > 0 and len(alt_trim) > 0:
        pos1 = deleted_start1
        ref = ref_trim
        alt = alt_trim

        #verify reference against FASTA
        ref_fa = fasta.fetch_seq_1based_closed(chrom_fasta, pos1, pos1 + len(ref) - 1).upper()
        if ref_fa != ref:
            raise ValueError(
                f"Reference mismatch for substitution/delins at {chrom_fasta}:{pos1}: "
                f"SPDI_REF={ref} FASTA_REF={ref_fa}"
            )
        return {
            "CHROM": chrom_fasta,
            "POS": pos1,
            "REF": ref,
            "ALT": alt
        }

    #pure insertion
    if len(ref_trim) == 0 and len(alt_trim) > 0:
        # insertion occurs between bases pos0 and pos0+1 in 1-based coordinates
        anchor_pos1 = pos0
        if anchor_pos1 < 1:
            raise ValueError("Insertion at very first base cannot be anchored safely for VCF.")
        anchor = fasta.fetch_base_1based(chrom_fasta, anchor_pos1).upper()
        return {
            "CHROM": chrom_fasta,
            "POS": anchor_pos1,
            "REF": anchor,
            "ALT": anchor + alt_trim
        }

    #pure deletion
    if len(ref_trim) > 0 and len(alt_trim) == 0:
        # deletion starts at deleted_start1, VCF needs anchor base before it
        anchor_pos1 = deleted_start1 - 1
        if anchor_pos1 < 1:
            raise ValueError("Deletion at first base cannot be anchored safely for VCF.")
        anchor = fasta.fetch_base_1based(chrom_fasta, anchor_pos1).upper()

        ref_deleted = fasta.fetch_seq_1based_closed(
            chrom_fasta,
            deleted_start1,
            deleted_start1 + len(ref_trim) - 1
        ).upper()

        if ref_deleted != ref_trim:
            raise ValueError(
                f"Reference mismatch for deletion at {chrom_fasta}:{deleted_start1}: "
                f"SPDI_REF={ref_trim} FASTA_REF={ref_deleted}"
            )

        return {
            "CHROM": chrom_fasta,
            "POS": anchor_pos1,
            "REF": anchor + ref_trim,
            "ALT": anchor
        }

    raise ValueError("Empty ref and alt after trimming, unsupported.")


def spliceai_supports_vcf_alleles(ref, alt):
    """
    SpliceAI README:
    only SNVs and simple INDELs where REF or ALT is a single base.
    """
    return (len(ref) == 1 or len(alt) == 1)


#Main
df = pd.read_csv(INPUT_CSV)
fasta = FastaReader(REFERENCE_FASTA)

included = []
excluded = []

for idx, row in df.iterrows():
    name = row.get("Name")
    variation_id = row.get("VariationID")
    variant_type = row.get("Variant type")
    chrom_raw = row.get("GRCh38Chromosome")
    chrom = chrom_to_fasta_name(chrom_raw)
    spdi = row.get("Canonical SPDI")

    if pd.isna(spdi):
        excluded.append({
            "row_index": idx,
            "VariationID": variation_id,
            "Name": name,
            "Variant type": variant_type,
            "Canonical SPDI": spdi,
            "reason": "missing_Canonical_SPDI"
        })
        continue

    if pd.isna(chrom):
        excluded.append({
            "row_index": idx,
            "VariationID": variation_id,
            "Name": name,
            "Variant type": variant_type,
            "Canonical SPDI": spdi,
            "reason": "missing_GRCh38Chromosome"
        })
        continue

    try:
        spdi_obj = parse_spdi(spdi)
        if spdi_obj is None:
            raise ValueError("bad_SPDI_format")

        v = spdi_to_vcf_record(spdi_obj, chrom, fasta)

        if not spliceai_supports_vcf_alleles(v["REF"], v["ALT"]):
            excluded.append({
                "row_index": idx,
                "VariationID": variation_id,
                "Name": name,
                "Variant type": variant_type,
                "Canonical SPDI": spdi,
                "reason": f"not_simple_for_spliceai_REFlen{len(v['REF'])}_ALTlen{len(v['ALT'])}"
            })
            continue

        ref_check = fasta.fetch_seq_1based_closed(
            v["CHROM"], v["POS"], v["POS"] + len(v["REF"]) - 1
        ).upper()
        if ref_check != v["REF"]:
            raise ValueError(
                f"final_VCF_REF_mismatch: FASTA={ref_check}, VCF_REF={v['REF']}"
            )

        variant_id = f"ClinVar_{variation_id}" if pd.notna(variation_id) else f"ClinVarRow_{idx}"

        included.append({
            "CHROM": v["CHROM"],
            "POS": int(v["POS"]),
            "ID": variant_id,
            "REF": v["REF"],
            "ALT": v["ALT"],
            "QUAL": ".",
            "FILTER": "PASS",
            "INFO": f"VariationID={variation_id};Name={str(name).replace(';', ',')};VariantType={str(variant_type).replace(';', ',')}",
            "VariationID": variation_id,
            "Name": name,
            "Variant type": variant_type,
            "Canonical SPDI": spdi
        })

    except Exception as e:
        excluded.append({
            "row_index": idx,
            "VariationID": variation_id,
            "Name": name,
            "Variant type": variant_type,
            "Canonical SPDI": spdi,
            "reason": str(e)
        })

included_df = pd.DataFrame(included)
excluded_df = pd.DataFrame(excluded)

if not included_df.empty:
    included_df = included_df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"]).copy()
    included_df = included_df.sort_values(["CHROM", "POS", "REF", "ALT"]).reset_index(drop=True)

#write VCF
#read fai for contig header
fai_file = REFERENCE_FASTA + ".fai"
contigs = []
with open(fai_file, "r") as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        contig_id = fields[0]
        contig_len = fields[1]
        contigs.append((contig_id, contig_len))

with open(OUT_VCF, "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    out.write(f"##reference={REFERENCE_FASTA}\n")
    for cid, clen in contigs:
        out.write(f"##contig=<ID={cid},length={clen}>\n")
    out.write('##INFO=<ID=VariationID,Number=1,Type=String,Description="ClinVar VariationID">\n')
    out.write('##INFO=<ID=Name,Number=1,Type=String,Description="ClinVar variant name">\n')
    out.write('##INFO=<ID=VariantType,Number=1,Type=String,Description="ClinVar variant type">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    if not included_df.empty:
        for _, r in included_df.iterrows():
            out.write(
                f"{r['CHROM']}\t{r['POS']}\t{r['ID']}\t{r['REF']}\t{r['ALT']}\t"
                f"{r['QUAL']}\t{r['FILTER']}\t{r['INFO']}\n"
            )

    if not included_df.empty:
        for _, r in included_df.iterrows():
            out.write(
                f"{r['CHROM']}\t{r['POS']}\t{r['ID']}\t{r['REF']}\t{r['ALT']}\t"
                f"{r['QUAL']}\t{r['FILTER']}\t{r['INFO']}\n"
            )

included_df.to_csv(OUT_INCLUDED, sep="\t", index=False)
excluded_df.to_csv(OUT_EXCLUDED, sep="\t", index=False)

with open(OUT_REPORT, "w") as out:
    out.write("ClinVar ABCA4 -> SpliceAI input VCF report\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"Input rows: {df.shape[0]}\n")
    out.write(f"Included VCF rows: {included_df.shape[0]}\n")
    out.write(f"Excluded rows: {excluded_df.shape[0]}\n\n")

    out.write("Included variant types:\n")
    if not included_df.empty:
        for k, v in included_df["Variant type"].value_counts().items():
            out.write(f"{k}\t{v}\n")
    out.write("\nExcluded reasons:\n")
    if not excluded_df.empty:
        for k, v in excluded_df["reason"].value_counts().items():
            out.write(f"{k}\t{v}\n")

print("Done.")
print("VCF:", OUT_VCF)
print("Included:", OUT_INCLUDED)
print("Excluded:", OUT_EXCLUDED)
print("Report:", OUT_REPORT)

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################