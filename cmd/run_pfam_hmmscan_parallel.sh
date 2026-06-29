########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

#!/usr/bin/env bash
set -euo pipefail

PROT_FA="./RetinaMerge.filtered_lite.orfanage.protein.fa"
PFAM_HMM="./pfamdb/Pfam-A.hmm"
OUTDIR="./pfam_analysis_parallel"

#16-thread machine recommended settings
N_JOBS=4
CPU_PER_JOB=4

LOGFILE="${OUTDIR}/run_pfam_hmmscan_parallel.log"
CHUNK_DIR="${OUTDIR}/chunks"
RESULT_DIR="${OUTDIR}/chunk_results"
MERGED_DOMTBLOUT="${OUTDIR}/pfam.domtblout"
MERGED_TBLOUT="${OUTDIR}/pfam.tblout"
SUMMARY_FILE="${OUTDIR}/pfam_run_summary.txt"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

log() {
  echo "[$(timestamp)] $*" | tee -a "${LOGFILE}"
}

die() {
  echo "[$(timestamp)] [ERROR] $*" | tee -a "${LOGFILE}" >&2
  exit 1
}

mkdir -p "${OUTDIR}" "${CHUNK_DIR}" "${RESULT_DIR}"
: > "${LOGFILE}"

log "===== Starting parallel Pfam hmmscan pipeline ====="
log "Working directory: $(pwd)"

command -v hmmscan >/dev/null 2>&1 || die "hmmscan not found in PATH"
command -v hmmpress >/dev/null 2>&1 || die "hmmpress not found in PATH"
command -v awk >/dev/null 2>&1 || die "awk not found"
command -v find >/dev/null 2>&1 || die "find not found"
command -v xargs >/dev/null 2>&1 || die "xargs not found"

[[ -f "${PROT_FA}" ]] || die "Protein FASTA not found: ${PROT_FA}"
[[ -s "${PROT_FA}" ]] || die "Protein FASTA is empty: ${PROT_FA}"
[[ -f "${PFAM_HMM}" ]] || die "Pfam HMM file not found: ${PFAM_HMM}"
[[ -s "${PFAM_HMM}" ]] || die "Pfam HMM file is empty: ${PFAM_HMM}"

SEQ_COUNT=$(grep -c '^>' "${PROT_FA}" || true)
[[ "${SEQ_COUNT}" -gt 0 ]] || die "No FASTA records detected"

log "Protein FASTA: ${PROT_FA}"
log "Pfam HMM: ${PFAM_HMM}"
log "Number of protein sequences: ${SEQ_COUNT}"
log "N_JOBS=${N_JOBS}"
log "CPU_PER_JOB=${CPU_PER_JOB}"

H3M="${PFAM_HMM}.h3m"
H3I="${PFAM_HMM}.h3i"
H3F="${PFAM_HMM}.h3f"
H3P="${PFAM_HMM}.h3p"

if [[ -f "${H3M}" && -f "${H3I}" && -f "${H3F}" && -f "${H3P}" ]]; then
  log "Pfam HMM database already pressed."
else
  log "Pressed database files not found. Running hmmpress ..."
  hmmpress "${PFAM_HMM}" >> "${LOGFILE}" 2>&1 || die "hmmpress failed"
  log "hmmpress completed."
fi

log "Splitting FASTA into ${N_JOBS} chunks ..."

rm -f "${CHUNK_DIR}"/chunk_*.fa

awk -v n="${N_JOBS}" -v outdir="${CHUNK_DIR}" '
BEGIN {
  seq = 0
}
(/^>/) {
  file_idx = seq % n
  outfile = sprintf("%s/chunk_%03d.fa", outdir, file_idx)
  seq++
}
{
  print >> outfile
}
' "${PROT_FA}"

CHUNK_COUNT=$(find "${CHUNK_DIR}" -maxdepth 1 -name 'chunk_*.fa' | wc -l)
log "Created ${CHUNK_COUNT} chunk FASTA files."

run_one_chunk() {
  local chunk_fa="$1"
  local base
  base=$(basename "${chunk_fa}" .fa)

  hmmscan \
    --cpu "${CPU_PER_JOB}" \
    --domtblout "${RESULT_DIR}/${base}.domtblout" \
    --tblout "${RESULT_DIR}/${base}.tblout" \
    -o "${RESULT_DIR}/${base}.out" \
    "${PFAM_HMM}" \
    "${chunk_fa}"
}

export PFAM_HMM RESULT_DIR CPU_PER_JOB
export -f run_one_chunk

log "Running hmmscan in parallel ..."
find "${CHUNK_DIR}" -maxdepth 1 -name 'chunk_*.fa' | sort | \
  xargs -I{} -P "${N_JOBS}" bash -c 'run_one_chunk "$@"' _ {}

log "All chunk hmmscan jobs finished."

log "Merging outputs ..."

first_dom=$(find "${RESULT_DIR}" -maxdepth 1 -name 'chunk_*.domtblout' | sort | head -n 1)
first_tbl=$(find "${RESULT_DIR}" -maxdepth 1 -name 'chunk_*.tblout' | sort | head -n 1)

[[ -n "${first_dom}" ]] || die "No chunk domtblout files found"
[[ -n "${first_tbl}" ]] || die "No chunk tblout files found"

cp "${first_dom}" "${MERGED_DOMTBLOUT}"
cp "${first_tbl}" "${MERGED_TBLOUT}"

for f in $(find "${RESULT_DIR}" -maxdepth 1 -name 'chunk_*.domtblout' | sort | tail -n +2); do
  grep -v '^#' "${f}" >> "${MERGED_DOMTBLOUT}"
done

for f in $(find "${RESULT_DIR}" -maxdepth 1 -name 'chunk_*.tblout' | sort | tail -n +2); do
  grep -v '^#' "${f}" >> "${MERGED_TBLOUT}"
done

DOMAIN_HIT_LINES=$(grep -vc '^#' "${MERGED_DOMTBLOUT}" || true)
SEQ_HIT_LINES=$(grep -vc '^#' "${MERGED_TBLOUT}" || true)

{
  echo "Pfam parallel hmmscan run summary"
  echo "================================="
  echo "Run time: $(timestamp)"
  echo "Working directory: $(pwd)"
  echo "Protein FASTA: ${PROT_FA}"
  echo "Pfam HMM: ${PFAM_HMM}"
  echo "Number of protein sequences: ${SEQ_COUNT}"
  echo "Number of chunks: ${CHUNK_COUNT}"
  echo "Parallel jobs: ${N_JOBS}"
  echo "CPU per job: ${CPU_PER_JOB}"
  echo "Sequence-level hit lines in merged tblout: ${SEQ_HIT_LINES}"
  echo "Domain-level hit lines in merged domtblout: ${DOMAIN_HIT_LINES}"
} > "${SUMMARY_FILE}"

log "Summary written to: ${SUMMARY_FILE}"
log "===== Parallel Pfam hmmscan pipeline finished successfully ====="

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################