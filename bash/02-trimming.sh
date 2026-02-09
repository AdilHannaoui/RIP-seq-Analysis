#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Trimming Module (parallel + pigz)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

# --------------------------
# Setup directories
# --------------------------
mkdir -p "$OUTPUT_DIR/fastqc_trimmed" "$OUTPUT_DIR/logs" "$FASTQ_TRIM"
cd "$WORKDIR"

# --------------------------
# Validate dependencies
# --------------------------
for cmd in parallel fastqc; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "ERROR: Required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# Java solo si vamos a usar JAR
if ! command -v trimmomatic &> /dev/null && ! command -v java &> /dev/null; then
    echo "ERROR: Neither 'trimmomatic' wrapper nor 'java' found in PATH" >&2
    exit 1
fi

# --------------------------
# Validate Trimmomatic (conda wrapper or JAR)
# --------------------------
if command -v trimmomatic &> /dev/null; then
    # Usar el wrapper de conda (más simple)
    TRIMMO_CMD="trimmomatic"
    ADAPTER_DIR="${CONDA_PREFIX}/share/trimmomatic/adapters"
elif [[ -f "${TRIMMO_JAR:-}" ]]; then
    # Fallback al JAR tradicional
    TRIMMO_CMD="java -jar $TRIMMO_JAR"
    ADAPTER_DIR="$(dirname "$TRIMMO_JAR")/../adapters"
else
    echo "ERROR: Trimmomatic not found (neither wrapper nor JAR)" >&2
    echo "  - Install via conda: conda install -c bioconda trimmomatic" >&2
    echo "  - Or set TRIMMO_JAR in config.sh" >&2
    exit 1
fi

ADAPTER_FILE="$ADAPTER_DIR/TruSeq3-SE.fa"
if [[ ! -f "$ADAPTER_FILE" ]]; then
    echo "ERROR: Adapter file not found at $ADAPTER_FILE" >&2
    echo "Available adapters in $ADAPTER_DIR:" >&2
    ls -1 "$ADAPTER_DIR"/*.fa 2>/dev/null || echo "  (none found)" >&2
    exit 1
fi

echo "Using Trimmomatic: $TRIMMO_CMD"
echo "Using adapters: $ADAPTER_FILE"

# --------------------------
# Detect FASTQ files (any extension)
# --------------------------
shopt -s nullglob
FASTQ_FILES=(
    "$FASTQ_DIR"/*.fastq 
    "$FASTQ_DIR"/*.fq
    "$FASTQ_DIR"/*.fastq.gz 
    "$FASTQ_DIR"/*.fq.gz
)
shopt -u nullglob

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in $FASTQ_DIR" >&2
    echo "Searched for: *.fastq, *.fq, *.fastq.gz, *.fq.gz" >&2
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} FASTQ files."

# --------------------------
# Decompress .gz files if pigz is available
# --------------------------
shopt -s nullglob
GZ_FILES=("$FASTQ_DIR"/*.fastq.gz "$FASTQ_DIR"/*.fq.gz)
shopt -u nullglob

if [[ ${#GZ_FILES[@]} -gt 0 ]]; then
    if command -v pigz &> /dev/null; then
        echo "Decompressing ${#GZ_FILES[@]} .gz files using pigz with $THREADS threads..."
        printf '%s\n' "${GZ_FILES[@]}" | parallel -j "$THREADS" pigz -d -p 1 {}
    else
        echo "WARNING: pigz not found, using gzip (slower)..." >&2
        for gz_file in "${GZ_FILES[@]}"; do
            echo "Decompressing $(basename "$gz_file")..."
            gzip -d "$gz_file"
        done
    fi
    
    # Refresh list after decompression
    shopt -s nullglob
    FASTQ_FILES=(
        "$FASTQ_DIR"/*.fastq 
        "$FASTQ_DIR"/*.fq
    )
    shopt -u nullglob
    
    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found after decompression" >&2
        exit 1
    fi
    
    echo "After decompression: ${#FASTQ_FILES[@]} FASTQ files ready."
fi

# --------------------------
# Function to trim a single FASTQ file
# --------------------------
run_trimming() {
    local FASTQ_FILE="$1"
    local SAMPLE_NAME EXTENSION
    
    # Extract sample name (remove .fastq, .fq, or any extension)
    SAMPLE_NAME=$(basename "$FASTQ_FILE")
    SAMPLE_NAME="${SAMPLE_NAME%.fastq}"
    SAMPLE_NAME="${SAMPLE_NAME%.fq}"
    
    local TRIMMED_FASTQ_FILE="$FASTQ_TRIM/${SAMPLE_NAME}_trimmed.fastq"
    local LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}.trimmomatic.log"
    
    echo ">>> Trimming $SAMPLE_NAME"
    
    if $TRIMMO_CMD SE -threads 1 \
        "$FASTQ_FILE" "$TRIMMED_FASTQ_FILE" \
        ILLUMINACLIP:"$ADAPTER_FILE:2:30:10" \
        SLIDINGWINDOW:4:20 MINLEN:20 -phred33 \
        > "$LOGFILE" 2>&1; then
        
        # Verify output file exists and is not empty
        if [[ ! -s "$TRIMMED_FASTQ_FILE" ]]; then
            echo "ERROR: Trimmomatic produced empty output for $SAMPLE_NAME" >&2
            return 1
        fi
        
        echo ">>> Finished trimming $SAMPLE_NAME (OK)"
    else
        echo "ERROR: Trimmomatic failed for $SAMPLE_NAME (check $LOGFILE)" >&2
        return 1
    fi
}

export -f run_trimming
export OUTPUT_DIR FASTQ_TRIM TRIMMO_CMD ADAPTER_FILE

# --------------------------
# Run Trimmomatic in parallel
# --------------------------
echo "Running Trimmomatic in parallel using $THREADS threads..."

if ! parallel -j "$THREADS" --halt soon,fail=1 run_trimming ::: "${FASTQ_FILES[@]}"; then
    echo "ERROR: Some Trimmomatic jobs failed. Check logs in $OUTPUT_DIR/logs" >&2
    exit 1
fi

# --------------------------
# Run FastQC on trimmed files
# --------------------------
shopt -s nullglob
TRIMMED_FILES=("$FASTQ_TRIM"/*_trimmed.fastq)
shopt -u nullglob

if [[ ${#TRIMMED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No trimmed files found in $FASTQ_TRIM" >&2
    exit 1
fi

echo "Running FastQC on ${#TRIMMED_FILES[@]} trimmed files..."

run_fastqc_trimmed() {
    local FASTQ_FILE="$1"
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)
    
    local LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}.fastqc.log"
    
    if fastqc "$FASTQ_FILE" \
        --threads 1 \
        --outdir "$OUTPUT_DIR/fastqc_trimmed" \
        > "$LOGFILE" 2>&1; then
        echo ">>> FastQC finished for $SAMPLE_NAME (OK)"
    else
        echo "ERROR: FastQC failed for $SAMPLE_NAME" >&2
        return 1
    fi
}

export -f run_fastqc_trimmed

if ! parallel -j "$THREADS" --halt soon,fail=1 run_fastqc_trimmed ::: "${TRIMMED_FILES[@]}"; then
    echo "ERROR: Some FastQC jobs failed. Check logs in $OUTPUT_DIR/logs" >&2
    exit 1
fi

echo "All trimming and FastQC analyses completed successfully."
echo "Trimmed FASTQ saved in: $FASTQ_TRIM"
echo "FastQC results saved in: $OUTPUT_DIR/fastqc_trimmed"
echo "Logs saved in: $OUTPUT_DIR/logs"