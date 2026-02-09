#!/usr/bin/env bash
set -euo pipefail

# ==========================
# featureCounts Quantification Module (parallel)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

# --------------------------
# Setup directories
# --------------------------
mkdir -p "$OUTPUT_DIR/featurecounts" "$OUTPUT_DIR/logs"
cd "$WORKDIR"

# --------------------------
# Validate dependencies
# --------------------------
if ! command -v featureCounts &> /dev/null; then
    echo "ERROR: featureCounts not found in PATH" >&2
    echo "Install via: conda install -c bioconda subread" >&2
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "ERROR: parallel not found in PATH" >&2
    exit 1
fi

# --------------------------
# Validate required paths
# --------------------------
if [[ -z "${GTF:-}" ]]; then
    echo "ERROR: GTF not defined in config.sh" >&2
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF file not found at: $GTF" >&2
    exit 1
fi

if [[ -z "${BAM_DIR:-}" ]]; then
    echo "ERROR: BAM_DIR not defined in config.sh" >&2
    exit 1
fi

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found at: $BAM_DIR" >&2
    exit 1
fi

echo "Using GTF: $GTF"
echo "Using BAM directory: $BAM_DIR"

# --------------------------
# Detect BAM files directly (more reliable than inferring from FASTQ)
# --------------------------
shopt -s nullglob
BAM_FILES=("$BAM_DIR"/*.bam)
shopt -u nullglob

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in $BAM_DIR" >&2
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM files to quantify."

# --------------------------
# Validate BAM files are indexed
# --------------------------
echo "Validating BAM files and indices..."
MISSING_INDEX=0

for bam in "${BAM_FILES[@]}"; do
    if [[ ! -f "${bam}.bai" ]]; then
        echo "WARNING: Index missing for $(basename "$bam"), creating it..." >&2
        if ! samtools index "$bam"; then
            echo "ERROR: Failed to index $(basename "$bam")" >&2
            ((MISSING_INDEX++))
        fi
    fi
done

if [[ $MISSING_INDEX -gt 0 ]]; then
    echo "ERROR: Failed to index $MISSING_INDEX BAM file(s)" >&2
    exit 1
fi

# --------------------------
# Function to run featureCounts
# --------------------------
run_featurecounts() {
    local BAM_FILE="$1"
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    local OUTFILE="$OUTPUT_DIR/featurecounts/${SAMPLE_NAME}_featurecounts.txt"
    local SUMMARY_FILE="${OUTFILE}.summary"
    local COUNTS_FILE="$OUTPUT_DIR/featurecounts/Counts_${SAMPLE_NAME}.txt"
    local LOGFILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}_featurecounts.log"
    
    echo ">>> Quantifying $SAMPLE_NAME"
    
    # Verify BAM file is not empty
    if [[ ! -s "$BAM_FILE" ]]; then
        echo "ERROR: BAM file is empty for $SAMPLE_NAME" >&2
        return 1
    fi
    
    # Run featureCounts
    if featureCounts \
        -T 1 \
        -s 2 \
        -a "$GTF" \
        -o "$OUTFILE" \
        "$BAM_FILE" \
        > "$LOGFILE" 2>&1; then
        
        # Verify output was created and is not empty
        if [[ ! -s "$OUTFILE" ]]; then
            echo "ERROR: featureCounts produced empty output for $SAMPLE_NAME" >&2
            return 1
        fi
        
        # Extract gene ID + counts (skip header with metadata lines)
        if awk 'NR > 1 && !/^#/ {print $1"\t"$7}' "$OUTFILE" > "$COUNTS_FILE"; then
            
            if [[ ! -s "$COUNTS_FILE" ]]; then
                echo "WARNING: No counts extracted for $SAMPLE_NAME (possible low mapping rate)" >&2
            fi
            
            echo ">>> Finished $SAMPLE_NAME (OK)"
        else
            echo "ERROR: Failed to extract counts for $SAMPLE_NAME" >&2
            return 1
        fi
        
    else
        echo "ERROR: featureCounts failed for $SAMPLE_NAME (check $LOGFILE)" >&2
        return 1
    fi
}

export -f run_featurecounts
export BAM_DIR OUTPUT_DIR GTF

# --------------------------
# Run featureCounts in parallel
# --------------------------
echo "Running featureCounts in parallel using $THREADS threads..."

if ! parallel -j "$THREADS" --halt soon,fail=1 run_featurecounts ::: "${BAM_FILES[@]}"; then
    echo "ERROR: Some featureCounts jobs failed. Check logs in $OUTPUT_DIR/logs" >&2
    exit 1
fi

# --------------------------
# Summary statistics
# --------------------------
echo ""
echo "=== featureCounts Summary ==="

shopt -s nullglob
COUNT_FILES=("$OUTPUT_DIR/featurecounts"/Counts_*.txt)
SUMMARY_FILES=("$OUTPUT_DIR/featurecounts"/*.summary)
shopt -u nullglob

echo "Count files generated: ${#COUNT_FILES[@]}"
echo ""

# Extract assignment statistics from summary files
if [[ ${#SUMMARY_FILES[@]} -gt 0 ]]; then
    echo "=== Assignment Statistics ==="
    for summary in "${SUMMARY_FILES[@]}"; do
        sample=$(basename "$summary" _featurecounts.txt.summary)
        if [[ -f "$summary" ]]; then
            assigned=$(awk 'NR==2 {print $2}' "$summary" 2>/dev/null || echo "0")
            total=$(awk 'NR>1 {sum+=$2} END {print sum}' "$summary" 2>/dev/null || echo "1")
            
            if [[ $total -gt 0 ]]; then
                percent=$(awk -v a="$assigned" -v t="$total" 'BEGIN {printf "%.2f", (a/t)*100}')
                echo "$sample: $assigned / $total assigned ($percent%)"
            fi
        fi
    done
fi

echo ""
echo "All featureCounts quantifications completed successfully."
echo "Full results saved in: $OUTPUT_DIR/featurecounts"
echo "Count matrices saved as: Counts_*.txt"
echo "Logs saved in: $OUTPUT_DIR/logs"