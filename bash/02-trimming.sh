#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Trimming Module (parallel + pigz)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/fastqc_trimmed"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$FASTQ_TRIM"

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
# Function to trim a single FASTQ file
# --------------------------
run_trimming() {
    FASTQ_FILE="$1"

    SAMPLE_NAME=$(basename "$FASTQ_FILE")
    SAMPLE_NAME="${SAMPLE_NAME%%.*}"   # remove .fastq or .fastq.gz

    TRIMMED_FASTQ_FILE="$FASTQ_TRIM/${SAMPLE_NAME}_trimmed.fastq"
    LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}.trimmomatic.log"

    echo ">>> Trimming $SAMPLE_NAME"

    java -jar "$TRIMMO_JAR" SE -threads 1 \
        "$FASTQ_FILE" "$TRIMMED_FASTQ_FILE" \
        ILLUMINACLIP:"$WORKDIR/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10" \
        SLIDINGWINDOW:4:20 MINLEN:20 -phred33 \
        > "$LOGFILE" 2>&1

    if [[ ! -s "$TRIMMED_FASTQ_FILE" ]]; then
        echo "ERROR: Trimmomatic failed for $SAMPLE_NAME"
        exit 1
    fi

    echo ">>> Finished trimming $SAMPLE_NAME"
}

export -f run_trimming
export OUTPUT_DIR FASTQ_TRIM TRIMMO_JAR WORKDIR

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
# Run Trimmomatic in parallel
# --------------------------
echo "Running Trimmomatic in parallel using $THREADS threads..."

parallel -j "$THREADS" run_trimming ::: "${FASTQ_FILES[@]}"

# --------------------------
# Run FastQC on trimmed files
# --------------------------
echo "Running FastQC on trimmed files..."

parallel -j "$THREADS" fastqc {} --outdir "$OUTPUT_DIR/fastqc_trimmed" \
    ">" "$OUTPUT_DIR/logs/{/.}_trimmed.fastqc.log" 2>&1 ::: "$FASTQ_TRIM"/*.fastq

echo "All trimming and FastQC analyses completed."
echo "Trimmed FASTQ saved in: $FASTQ_TRIM"
echo "FastQC results saved in: $OUTPUT_DIR/fastqc_trimmed"
echo "Logs saved in: $OUTPUT_DIR/logs"
