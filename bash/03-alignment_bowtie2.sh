#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Bowtie2 Alignment Module (parallel + pigz)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

# --------------------------
# Setup directories
# --------------------------
mkdir -p "$OUTPUT_DIR/bowtie2" "$OUTPUT_DIR/visualization" "$OUTPUT_DIR/logs"
cd "$WORKDIR"

# --------------------------
# Validate dependencies
# --------------------------
for cmd in bowtie2 samtools bedtools parallel; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "ERROR: Required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# --------------------------
# Validate Bowtie2 index
# --------------------------
if [[ -z "${BOWTIE2_INDEX:-}" ]]; then
    echo "ERROR: BOWTIE2_INDEX not defined in config.sh" >&2
    exit 1
fi

# Check if index files exist (at least .1.bt2 file)
if ! ls "${BOWTIE2_INDEX}".*.bt2 &> /dev/null && ! ls "${BOWTIE2_INDEX}".*.bt2l &> /dev/null; then
    echo "ERROR: Bowtie2 index not found at: $BOWTIE2_INDEX" >&2
    echo "Expected files like: ${BOWTIE2_INDEX}.1.bt2 or ${BOWTIE2_INDEX}.1.bt2l" >&2
    exit 1
fi

echo "Using Bowtie2 index: $BOWTIE2_INDEX"

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
# Decompress .gz files if needed
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
# Function to align a single FASTQ file
# --------------------------
run_alignment() {
    local FASTQ_FILE="$1"
    local SAMPLE_NAME
    
    # Extract sample name (remove extensions)
    SAMPLE_NAME=$(basename "$FASTQ_FILE")
    SAMPLE_NAME="${SAMPLE_NAME%.fastq}"
    SAMPLE_NAME="${SAMPLE_NAME%.fq}"
    
    local BAM_FILE="$OUTPUT_DIR/bowtie2/${SAMPLE_NAME}.bam"
    local BEDGRAPH_FILE="$OUTPUT_DIR/visualization/${SAMPLE_NAME}.bedgraph"
    local LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}.bowtie2.log"
    
    echo ">>> Aligning $SAMPLE_NAME"
    
    # Bowtie2 + Samtools sort (usando pipe para evitar archivos intermedios)
    if bowtie2 -p 1 --very-fast -x "$BOWTIE2_INDEX" -U "$FASTQ_FILE" 2> "$LOGFILE" \
        | samtools sort -@ 1 -o "$BAM_FILE" -; then
        
        # Verify BAM file was created and is not empty
        if [[ ! -s "$BAM_FILE" ]]; then
            echo "ERROR: Empty BAM file produced for $SAMPLE_NAME" >&2
            return 1
        fi
        
        # Index BAM file
        if ! samtools index "$BAM_FILE"; then
            echo "ERROR: Failed to index BAM for $SAMPLE_NAME" >&2
            return 1
        fi
        
        echo ">>> Generating BedGraph for $SAMPLE_NAME"
        
        # BedGraph generation (sorted)
        if bedtools genomecov -ibam "$BAM_FILE" -bg \
            | sort -k1,1 -k2,2n > "$BEDGRAPH_FILE"; then
            
            if [[ ! -s "$BEDGRAPH_FILE" ]]; then
                echo "WARNING: Empty BedGraph produced for $SAMPLE_NAME" >&2
            fi
        else
            echo "ERROR: Bedtools failed for $SAMPLE_NAME" >&2
            return 1
        fi
        
        echo ">>> Finished alignment for $SAMPLE_NAME (OK)"
    else
        echo "ERROR: Bowtie2 or Samtools failed for $SAMPLE_NAME (check $LOGFILE)" >&2
        return 1
    fi
}

export -f run_alignment
export OUTPUT_DIR BOWTIE2_INDEX

# --------------------------
# Run Bowtie2 alignment in parallel
# --------------------------
echo "Running Bowtie2 alignment in parallel using $THREADS threads..."

if ! parallel -j "$THREADS" --halt soon,fail=1 run_alignment ::: "${FASTQ_FILES[@]}"; then
    echo "ERROR: Some alignment jobs failed. Check logs in $OUTPUT_DIR/logs" >&2
    exit 1
fi

# --------------------------
# Summary statistics
# --------------------------
echo ""
echo "=== Alignment Summary ==="

shopt -s nullglob
BAM_FILES=("$OUTPUT_DIR/bowtie2"/*.bam)
BEDGRAPH_FILES=("$OUTPUT_DIR/visualization"/*.bedgraph)
shopt -u nullglob

echo "BAM files generated: ${#BAM_FILES[@]}"
echo "BedGraph files generated: ${#BEDGRAPH_FILES[@]}"

# Extract alignment rates from logs
echo ""
echo "=== Alignment Rates ==="
for log in "$OUTPUT_DIR/logs"/*.bowtie2.log; do
    if [[ -f "$log" ]]; then
        sample=$(basename "$log" .bowtie2.log)
        rate=$(grep "overall alignment rate" "$log" 2>/dev/null | head -n1 || echo "N/A")
        echo "$sample: $rate"
    fi
done

echo ""
echo "All alignments completed successfully."
echo "BAM files saved in: $OUTPUT_DIR/bowtie2"
echo "BedGraph files saved in: $OUTPUT_DIR/visualization"
echo "Logs saved in: $OUTPUT_DIR/logs"