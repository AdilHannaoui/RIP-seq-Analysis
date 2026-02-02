#!/usr/bin/env bash
set -euo pipefail

# ==========================
# featureCounts Quantification Module (parallel)
# Author: Adil Hannaoui Anaaoui
# ==========================
# Load global config
source "$(dirname "$0")/config.sh"

mkdir -p "$OUTPUT_DIR/featurecounts"
mkdir -p "$OUTPUT_DIR/logs"

cd "$WORKDIR"

# --------------------------
# Detect FASTQ files (to infer sample names)
# --------------------------
FASTQ_FILES=($(ls "$FASTQ_TRIM"/*.fastq 2>/dev/null || true))

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "No FASTQ files found in $FASTQ_TRIM"
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} samples."

# --------------------------
# Function to run featureCounts
# --------------------------
run_featurecounts() {
    FASTQ_FILE="$1"
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)

    BAM_FILE="$BAM_DIR/${SAMPLE_NAME}.bam"
    OUTFILE="$OUTPUT_DIR/featurecounts/${SAMPLE_NAME}_featurecounts.txt"
    LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}_featurecounts.log"

    echo ">>> Quantifying $SAMPLE_NAME"

    if [[ ! -f "$BAM_FILE" ]]; then
        echo "ERROR: BAM file not found for $SAMPLE_NAME" >&2
        exit 1
    fi

    # Run featureCounts
    featureCounts \
        -T 1 \
        -S 2 \
        -a "$GTF_FILE" \
        -o "$OUTFILE" \
        "$BAM_FILE" \
        > "$LOGFILE" 2>&1

    if [[ ! -s "$OUTFILE" ]]; then
        echo "ERROR: featureCounts failed for $SAMPLE_NAME" >&2
        exit 1
    fi

    # Extract gene ID + counts
    cut -f1,7 "$OUTFILE" | tail -n +2 > "$OUTPUT_DIR/featurecounts/Counts_${SAMPLE_NAME}.txt"

    echo ">>> Finished $SAMPLE_NAME"
}

export -f run_featurecounts
export BAM_DIR OUTPUT_DIR GTF_FILE

# --------------------------
# Run featureCounts in parallel
# --------------------------
echo "Running featureCounts in parallel using $THREADS threads..."

parallel -j "$THREADS" run_featurecounts ::: "${FASTQ_FILES[@]}"

echo "All featureCounts quantifications completed."
echo "Results saved in: $OUTPUT_DIR/featurecounts"
echo "Logs saved in: $OUTPUT_DIR/logs"
