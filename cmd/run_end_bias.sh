########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env bash
set -euo pipefail
#FLNC 5'/3' end bias analysis
#Input files
BAM=/mnt/d/bioinformatic/isoseq/processed/from_Ray/proActiv/RayMerged_bamfile/RetinaMerge.dedup.mapped.bam
ANN=/mnt/d/bioinformatic/isoseq/processed/from_Ray/GRCh38p14/GRCh38.p14.v44/gencode.v44.chr_patch_hapl_scaff.annotation.gff3

#Output directory
OUT=.

#Generate TSS and TES BED files from transcript models
#clean possible leftovers
rm -f "$OUT/tss.bed" "$OUT/tes.bed"

awk -F'\t' -v OUT="$OUT" '
BEGIN{OFS="\t"}

function get_attr(attr, key,   r, m){
  r = key "=([^;]+)"
  if (match(attr, r, m)) return m[1]
  return "NA"
}

# Accept both "transcript" and "mRNA" features (GENCODE GFF3)
$0 !~ /^#/ && ($3=="transcript" || $3=="mRNA"){
  chr=$1; start=$4; end=$5; strand=$7; attr=$9;

  # Transcript ID (GFF3: ID=)
  tid = get_attr(attr, "ID")

  # Gene ID (prefer gene_id=, fall back to Parent=)
  gid = get_attr(attr, "gene_id")
  if (gid=="NA") gid = get_attr(attr, "Parent")

  t_start0 = start - 1
  t_end0   = end

  if (strand=="+") {
    print chr, t_start0, t_start0+1, tid"|"gid, 0, strand >> (OUT "/tss.bed")
    print chr, t_end0-1, t_end0,     tid"|"gid, 0, strand >> (OUT "/tes.bed")
  } else if (strand=="-") {
    print chr, t_end0-1, t_end0,     tid"|"gid, 0, strand >> (OUT "/tss.bed")
    print chr, t_start0, t_start0+1, tid"|"gid, 0, strand >> (OUT "/tes.bed")
  }
}
' "$ANN"

#sort TSS / TES BEDs
sort -k1,1 -k2,2n "$OUT/tss.bed" > "$OUT/tss.sorted.bed"
sort -k1,1 -k2,2n "$OUT/tes.bed" > "$OUT/tes.sorted.bed"

#Sort and index BAM
samtools sort -o "$OUT/flnc.sorted.bam" "$BAM"
samtools index "$OUT/flnc.sorted.bam"

#Convert BAM to BED6 and sort reads
bedtools bamtobed -i "$OUT/flnc.sorted.bam" > "$OUT/reads.bed"
sort -k1,1 -k2,2n "$OUT/reads.bed" > "$OUT/reads.sorted.bed"

#Generate read-level 5' and 3' end BEDs 
#5' end: + strand = start; - strand = end-1
awk 'BEGIN{OFS="\t"}
{
  chr=$1; s=$2; e=$3; name=$4; strand=$6;
  if (strand=="+")      p=s;
  else if (strand=="-") p=e-1;
  else next;
  print chr, p, p+1, name, 0, strand;
}' "$OUT/reads.sorted.bed" > "$OUT/read5p.bed"

#3' end: + strand = end-1; - strand = start
awk 'BEGIN{OFS="\t"}
{
  chr=$1; s=$2; e=$3; name=$4; strand=$6;
  if (strand=="+")      p=e-1;
  else if (strand=="-") p=s;
  else next;
  print chr, p, p+1, name, 0, strand;
}' "$OUT/reads.sorted.bed" > "$OUT/read3p.bed"

sort -k1,1 -k2,2n "$OUT/read5p.bed" > "$OUT/read5p.sorted.bed"
sort -k1,1 -k2,2n "$OUT/read3p.bed" > "$OUT/read3p.sorted.bed"

#Distances to nearest annotated TSS / TES
bedtools closest -s -D a \
  -a "$OUT/read5p.sorted.bed" \
  -b "$OUT/tss.sorted.bed" \
  > "$OUT/read5p_vs_tss.closest.txt"

bedtools closest -s -D a \
  -a "$OUT/read3p.sorted.bed" \
  -b "$OUT/tes.sorted.bed" \
  > "$OUT/read3p_vs_tes.closest.txt"

#Generate gene body spans from GFF3
awk -F'\t' -v OUT="$OUT" '
BEGIN{OFS="\t"}

function get_attr(attr, key,   r, m){
  r = key "=([^;]+)"
  if (match(attr, r, m)) return m[1]
  return "NA"
}

$0 !~ /^#/ && $3=="gene"{
  chr=$1; start=$4; end=$5; strand=$7; attr=$9;

  # Prefer gene_id=..., otherwise use ID=...
  gid   = get_attr(attr, "gene_id")
  if (gid=="NA") gid = get_attr(attr, "ID")

  # Prefer gene_name=..., otherwise use Name=...
  gname = get_attr(attr, "gene_name")
  if (gname=="NA") gname = get_attr(attr, "Name")

  print chr, start-1, end, gid"|"gname, 0, strand
}
' "$ANN" > "$OUT/genes.bed"

sort -k1,1 -k2,2n "$OUT/genes.bed" > "$OUT/genes.sorted.bed"

#Assign each read to one gene
bedtools intersect -s -wa -wb \
  -a "$OUT/reads.sorted.bed" \
  -b "$OUT/genes.sorted.bed" \
  > "$OUT/read_gene.intersect.txt"

#Pick the gene with maximum overlap for each read.
#Output columns: read_name, chr, read_start, read_end, strand, gene_id, gene_name, gene_start, gene_end
awk 'BEGIN{OFS="\t"}
{
  rchr=$1; rs=$2; re=$3; rname=$4; rstrand=$6;
  gs=$8; ge=$9; gid=$10;

  os = (rs>gs?rs:gs);
  oe = (re<ge?re:ge);
  ov = oe-os;
  if(ov<=0) next;

  if(!(rname in best) || ov>best[rname]){
    best[rname]=ov;
    rec[rname]=rname"\t"rchr"\t"rs"\t"re"\t"rstrand"\t"gid"\t"gs"\t"ge;
  }
}
END{
  for(k in rec) print rec[k];
}' "$OUT/read_gene.intersect.txt" > "$OUT/read_best_gene.tsv"

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################