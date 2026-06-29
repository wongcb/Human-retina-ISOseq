########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env bash
set -euo pipefail

#paths
ROOT="/mnt/d/bioinformatic/isoseq/processed/from_Ray/original_from_Ray"

QUERY_GFF="${ROOT}/ray_pigeon_output/classify.retina.merge.GRCh38.p14.with.cage.polya.FL/RetinaMerge.sorted.filtered_lite.gff"

REF_GTF_GZ="/mnt/d/bioinformatic/isoseq/processed/from_Ray/GRCh38p14/GRCh38.p14.v44/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz"
REF_GENOME="/mnt/d/bioinformatic/isoseq/processed/from_Ray/GRCh38p14/GRCh38.p14.genome.fa"

OUTDIR="${ROOT}/orf_prediction_orfanage_filtered_lite"
TMPDIR="${OUTDIR}/tmp"
LOG="${OUTDIR}/run.log"

REF_GTF="${TMPDIR}/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
ORF_OUT="${OUTDIR}/RetinaMerge.filtered_lite.orfanage.gtf"
ORF_STATS="${OUTDIR}/RetinaMerge.filtered_lite.orfanage.stats.tsv"

QUERY_TRANSCRIPTS_FA="${OUTDIR}/RetinaMerge.filtered_lite.transcripts.fa"
ORF_CDS_FA="${OUTDIR}/RetinaMerge.filtered_lite.orfanage.cds.fa"
ORF_PROTEIN_FA="${OUTDIR}/RetinaMerge.filtered_lite.orfanage.protein.fa"

THREADS=12

mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"

{
  echo "=== ORFanage pipeline started ==="
  date
  echo "QUERY_GFF=${QUERY_GFF}"
  echo "REF_GTF_GZ=${REF_GTF_GZ}"
  echo "REF_GENOME=${REF_GENOME}"
  echo "OUTDIR=${OUTDIR}"
  echo "THREADS=${THREADS}"
} | tee "${LOG}"

#check dependencies
for exe in orfanage gffread gzip python3; do
  if ! command -v "${exe}" >/dev/null 2>&1; then
    echo "Required executable not found in PATH: ${exe}" | tee -a "${LOG}"
    exit 1
  fi
done

for f in "${QUERY_GFF}" "${REF_GTF_GZ}" "${REF_GENOME}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Missing required file: ${f}" | tee -a "${LOG}"
    exit 1
  fi
done

#decompress reference GTF
echo "Decompress reference GTF" | tee -a "${LOG}"

gzip -cd "${REF_GTF_GZ}" > "${REF_GTF}"

echo "Reference GTF ready: ${REF_GTF}" | tee -a "${LOG}"

#run ORFanage
echo "Run ORFanage" | tee -a "${LOG}"

orfanage \
  --query "${QUERY_GFF}" \
  --output "${ORF_OUT}" \
  --reference "${REF_GENOME}" \
  --cleanq \
  --mode BEST \
  --threads "${THREADS}" \
  --stats "${ORF_STATS}" \
  "${REF_GTF}" \
  2>&1 | tee -a "${LOG}"

echo "ORFanage annotation written: ${ORF_OUT}" | tee -a "${LOG}"
echo "ORFanage stats written: ${ORF_STATS}" | tee -a "${LOG}"

#export query transcript FASTA
echo "Export query transcript FASTA from filtered_lite GFF" | tee -a "${LOG}"

gffread \
  "${QUERY_GFF}" \
  -g "${REF_GENOME}" \
  -w "${QUERY_TRANSCRIPTS_FA}" \
  2>&1 | tee -a "${LOG}"

echo "Query transcript FASTA written: ${QUERY_TRANSCRIPTS_FA}" | tee -a "${LOG}"

#export CDS and protein FASTA from ORFanage output
echo "Export CDS and protein FASTA from ORFanage output" | tee -a "${LOG}"

gffread \
  "${ORF_OUT}" \
  -g "${REF_GENOME}" \
  -x "${ORF_CDS_FA}" \
  -y "${ORF_PROTEIN_FA}" \
  2>&1 | tee -a "${LOG}"

echo "CDS FASTA written: ${ORF_CDS_FA}" | tee -a "${LOG}"
echo "Protein FASTA written: ${ORF_PROTEIN_FA}" | tee -a "${LOG}"

#basic checks
echo "Basic checks" | tee -a "${LOG}"

python3 - << 'PY'
from pathlib import Path

query_fa = Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/original_from_Ray/orf_prediction_orfanage_filtered_lite/RetinaMerge.filtered_lite.transcripts.fa")
cds_fa   = Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/original_from_Ray/orf_prediction_orfanage_filtered_lite/RetinaMerge.filtered_lite.orfanage.cds.fa")
pep_fa   = Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/original_from_Ray/orf_prediction_orfanage_filtered_lite/RetinaMerge.filtered_lite.orfanage.protein.fa")
orf_gtf  = Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/original_from_Ray/orf_prediction_orfanage_filtered_lite/RetinaMerge.filtered_lite.orfanage.gtf")
stats_tsv= Path("/mnt/d/bioinformatic/isoseq/processed/from_Ray/original_from_Ray/orf_prediction_orfanage_filtered_lite/RetinaMerge.filtered_lite.orfanage.stats.tsv")

def count_fasta(path):
    if not path.exists():
        return None
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n

print("query_transcripts_fasta_records:", count_fasta(query_fa))
print("orf_cds_fasta_records:", count_fasta(cds_fa))
print("orf_protein_fasta_records:", count_fasta(pep_fa))
print("orf_gtf_exists:", orf_gtf.exists())
print("stats_exists:", stats_tsv.exists())
PY

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################