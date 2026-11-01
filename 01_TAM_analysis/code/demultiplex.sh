#!/usr/bin/env bash
set -euo pipefail
# Supply the input read, output directory, and (optionally) executable/threads.
: "${R1_FASTQ:?Set R1_FASTQ to the pooled R1 FASTQ}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to the demultiplexing output directory}"
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BARCODES=${BARCODES:-"$SCRIPT_DIR/../metadata/barcodes.csv"}
ULTRAPLEX=${ULTRAPLEX:-ultraplex}
THREADS=${THREADS:-8}
"$ULTRAPLEX" -i "$R1_FASTQ" -b "$BARCODES" -d "$OUTPUT_DIR" -t "$THREADS" -dbr
