#!/bin/bash
set -euo pipefail

DEMUX_DIR="${1:-analysis_minibwa_bam/demux}"
CLEAN_DIR="${2:-analysis_minibwa_bam/clean}"
THREADS="${3:-4}"
JOBS="${4:-4}"
FASTP_BIN="${FASTP_BIN:-fastp}"
FASTP_LIB_DIR="${FASTP_LIB_DIR:-}"
PCR_ENV_LIB_DIR="${PCR_ENV_LIB_DIR:-}"

mkdir -p "${CLEAN_DIR}"
[[ -z "$FASTP_LIB_DIR" ]] || export LD_LIBRARY_PATH="$FASTP_LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
[[ -z "$PCR_ENV_LIB_DIR" ]] || export LD_LIBRARY_PATH="$PCR_ENV_LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

samples=$(find "${DEMUX_DIR}" -maxdepth 1 \( -name '*_R1.fastq' -o -name '*_R1.fastq.gz' \) -printf '%f\n' | sed -E 's/_R1.fastq(.gz)?$//' | sort)

run_one() {
  local sample="$1"
  local in1="${DEMUX_DIR}/${sample}_R1.fastq"
  local in2="${DEMUX_DIR}/${sample}_R2.fastq"
  if [[ ! -f "${in1}" ]]; then
    in1="${in1}.gz"
    in2="${in2}.gz"
  fi
  if [[ -s "${CLEAN_DIR}/${sample}_R1.clean.fastq.gz" && -s "${CLEAN_DIR}/${sample}_R2.clean.fastq.gz" ]]; then
    echo "fastp skip ${sample}"
    return
  fi
  echo "fastp start ${sample}"
  "${FASTP_BIN}" \
    --in1 "${in1}" \
    --in2 "${in2}" \
    --out1 "${CLEAN_DIR}/${sample}_R1.clean.fastq.gz" \
    --out2 "${CLEAN_DIR}/${sample}_R2.clean.fastq.gz" \
    --thread "${THREADS}" \
    --json "${CLEAN_DIR}/${sample}.fastp.json" \
    --html "${CLEAN_DIR}/${sample}.fastp.html"
  echo "fastp done ${sample}"
}

export -f run_one
export DEMUX_DIR CLEAN_DIR THREADS FASTP_BIN LD_LIBRARY_PATH
printf '%s\n' ${samples} | xargs -n 1 -P "${JOBS}" bash -lc 'run_one "$@"' _
