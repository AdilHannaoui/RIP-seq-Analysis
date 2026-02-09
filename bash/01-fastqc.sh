#!/usr/bin/env bash
set -euo pipefail

# ==========================
# FastQC QC Module (parallel + pigz)
# Author: Adil Hannaoui Anaaoui
# ==========================
# Load global config
source "$(dirname "$0")/config.sh"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/fastqc_results"
mkdir -p "$OUTPUT_DIR/logs"

cd "$WORKDIR"

# --------------------------
# Detect FASTQ files (.fastq or .fastq.gz)
# --------------------------
FASTQ_FILES=($(ls "$FASTQ_DIR"/*.fastq "$FASTQ_DIR"/*.fastq.gz 2>/dev/null || true))

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} FASTQ files."

# --------------------------
# Function to run FastQC on a single file
# --------------------------
run_fastqc() {
    FASTQ_FILE="$1"
    SAMPLE_NAME=$(basename "$FASTQ_FILE")
    SAMPLE_NAME="${SAMPLE_NAME%%.*}"   # remove .fastq or .fastq.gz

    LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}.fastqc.log"

    echo ">>> Running FastQC for $SAMPLE_NAME"

    fastqc "$FASTQ_FILE" \
        --threads 1 \
        --outdir "$OUTPUT_DIR/fastqc_results" \
        > "$LOGFILE" 2>&1

    echo ">>> Finished $SAMPLE_NAME"
}

export -f run_fastqc
export OUTPUT_DIR

# --------------------------
# Optional: decompress gz files using pigz (parallel gzip)
# --------------------------
echo "Checking for compressed FASTQ files..."

for f in "${FASTQ_FILES[@]}"; do
    if [[ "$f" == *.gz ]]; then
        echo "Decompressing $f using pigz..."
        pigz -d -p "$THREADS" "$f"
    fi
done

# Refresh list after decompression
FASTQ_FILES=($(ls "$FASTQ_DIR"/*.fastq))

# --------------------------
# Run FastQC in parallel
# --------------------------
echo "Running FastQC in parallel using $THREADS threads..."

parallel -j "$THREADS" run_fastqc ::: "${FASTQ_FILES[@]}"

echo "All FastQC analyses completed."
echo "Results saved in: $OUTPUT_DIR/fastqc_results"
echo "Logs saved in: $OUTPUT_DIR/logs"
