#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Bowtie2 Alignment Module (parallel + pigz)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

mkdir -p "$OUTPUT_DIR/bowtie2"
mkdir -p "$OUTPUT_DIR/visualization"
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
# Function to align a single FASTQ file
# --------------------------
run_alignment() {
    FASTQ_FILE="$1"

    SAMPLE_NAME=$(basename "$FASTQ_FILE")
    SAMPLE_NAME="${SAMPLE_NAME%%.*}"   # remove .fastq or .fastq.gz

    BAM_FILE="$OUTPUT_DIR/bowtie2/${SAMPLE_NAME}_trimmed.bam"
    BEDGRAPH_FILE="$OUTPUT_DIR/visualization/${SAMPLE_NAME}.bedgraph"
    LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}.bowtie2.log"

    echo ">>> Aligning $SAMPLE_NAME"

    # Bowtie2 + Samtools sort
    if ! bowtie2 -p 1 -x "$BOWTIE2_INDEX" -U "$FASTQ_FILE" \
        2> "$LOGFILE" | samtools sort -@ 1 -o "$BAM_FILE"; then
        echo "ERROR: Bowtie2 or Samtools failed for $SAMPLE_NAME"
        exit 1
    fi

    samtools index "$BAM_FILE"

    # BedGraph generation
    if ! bedtools genomecov -ibam "$BAM_FILE" -bg > "$BEDGRAPH_FILE"; then
        echo "ERROR: Bedtools failed for $SAMPLE_NAME"
        exit 1
    fi

    echo ">>> Finished alignment for $SAMPLE_NAME"
}

export -f run_alignment
export OUTPUT_DIR BOWTIE2_INDEX

# --------------------------
# Optional: decompress gz files using pigz
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
# Run Bowtie2 alignment in parallel
# --------------------------
echo "Running Bowtie2 alignment in parallel using $THREADS threads..."

parallel -j "$THREADS" run_alignment ::: "${FASTQ_FILES[@]}"

echo "All alignments completed."
echo "BAM files saved in: $OUTPUT_DIR/bowtie2"
echo "BedGraph files saved in: $OUTPUT_DIR/visualization"
echo "Logs saved in: $OUTPUT_DIR/logs"
