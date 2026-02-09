#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq MACS3 Peak Calling (parallel)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

# --------------------------
# Setup directories
# --------------------------
mkdir -p "$OUTPUT_DIR/macs3" "$OUTPUT_DIR/logs"
cd "$WORKDIR"

# --------------------------
# Validate dependencies
# --------------------------
if ! command -v macs3 &> /dev/null; then
    echo "ERROR: macs3 not found in PATH" >&2
    echo "Install via: conda install -c bioconda macs3" >&2
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "ERROR: parallel not found in PATH" >&2
    exit 1
fi

# --------------------------
# Validate required paths and variables
# --------------------------
if [[ -z "${OUTPUT_DIR:-}" ]]; then
    echo "ERROR: OUTPUT_DIR not defined in config.sh" >&2
    exit 1
fi

if [[ -z "${MAX_REPLICAS:-}" ]]; then
    echo "ERROR: MAX_REPLICAS not defined in config.sh" >&2
    exit 1
fi

BAM_DIR="$OUTPUT_DIR/bowtie2"

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found at: $BAM_DIR" >&2
    exit 1
fi

echo "Using BAM directory: $BAM_DIR"
echo "Maximum replicates: $MAX_REPLICAS"

# --------------------------
# Function to run MACS3 for a single IP/IN pair
# --------------------------
run_macs3() {
    local IP_FILE="$1"
    local REP="$2"
    local BASENAME SAMPLE_NAME IN_FILE
    
    BASENAME=$(basename "$IP_FILE" .bam)
    
    # Extract sample name (remove _IP{rep} suffix, case-insensitive)
    SAMPLE_NAME=$(echo "$BASENAME" | sed -E "s/_[Ii][Pp]${REP}.*//")
    
    # Find matching IN file (case-insensitive)
    shopt -s nullglob nocaseglob
    local IN_CANDIDATES=("$BAM_DIR/${SAMPLE_NAME}_"[Ii][Nn]${REP}*.bam)
    shopt -u nullglob nocaseglob
    
    if [[ ${#IN_CANDIDATES[@]} -eq 0 ]]; then
        echo "ERROR: No IN file found for sample $SAMPLE_NAME (rep $REP)" >&2
        echo "  Expected pattern: ${BAM_DIR}/${SAMPLE_NAME}_IN${REP}*.bam" >&2
        return 1
    fi
    
    if [[ ${#IN_CANDIDATES[@]} -gt 1 ]]; then
        echo "WARNING: Multiple IN files found for $SAMPLE_NAME (rep $REP), using first match" >&2
    fi
    
    IN_FILE="${IN_CANDIDATES[0]}"
    
    # Verify both files exist and are not empty
    if [[ ! -s "$IP_FILE" ]]; then
        echo "ERROR: IP file is empty or missing: $IP_FILE" >&2
        return 1
    fi
    
    if [[ ! -s "$IN_FILE" ]]; then
        echo "ERROR: IN file is empty or missing: $IN_FILE" >&2
        return 1
    fi
    
    local OUTPUT_PREFIX="${SAMPLE_NAME}_rep${REP}"
    local LOGFILE="$OUTPUT_DIR/logs/macs3_${OUTPUT_PREFIX}.log"
    
    echo ">>> Running MACS3 for $SAMPLE_NAME (rep $REP)"
    echo "    IP: $(basename "$IP_FILE")"
    echo "    IN: $(basename "$IN_FILE")"
    
    if macs3 callpeak \
        -t "$IP_FILE" \
        -c "$IN_FILE" \
        --format BAM \
        --name "${OUTPUT_PREFIX}" \
        --outdir "$OUTPUT_DIR/macs3" \
        --pvalue 0.1 \
        --call-summits \
        --nomodel \
        --extsize 150 \
        --gsize 1.2e7 \
        --keep-dup all \
        > "$LOGFILE" 2>&1; then
        
        # Verify peak file was created
        local PEAK_FILE="$OUTPUT_DIR/macs3/${OUTPUT_PREFIX}_peaks.narrowPeak"
        
        if [[ -f "$PEAK_FILE" ]]; then
            local PEAK_COUNT
            PEAK_COUNT=$(wc -l < "$PEAK_FILE")
            echo ">>> Finished MACS3 for $SAMPLE_NAME (rep $REP): $PEAK_COUNT peaks called"
        else
            echo "WARNING: No peaks file generated for $SAMPLE_NAME (rep $REP)" >&2
        fi
    else
        echo "ERROR: MACS3 failed for $SAMPLE_NAME (rep $REP) - check $LOGFILE" >&2
        return 1
    fi
}

export -f run_macs3
export OUTPUT_DIR BAM_DIR

# --------------------------
# Build list of IP files for all replicates
# --------------------------
echo "Searching for IP BAM files..."

IP_LIST=()

for REP in $(seq 1 "$MAX_REPLICAS"); do
    shopt -s nullglob nocaseglob
    IP_FILES=("$BAM_DIR"/*_[Ii][Pp]${REP}*.bam)
    shopt -u nullglob nocaseglob
    
    for IP_FILE in "${IP_FILES[@]}"; do
        IP_LIST+=("$IP_FILE:$REP")
    done
done

if [[ ${#IP_LIST[@]} -eq 0 ]]; then
    echo "ERROR: No IP BAM files found in $BAM_DIR" >&2
    echo "Searched pattern: *_IP{1..$MAX_REPLICAS}*.bam (case-insensitive)" >&2
    exit 1
fi

echo "Found ${#IP_LIST[@]} IP files for MACS3 peak calling."

# Print IP-IN pairs that will be processed
echo ""
echo "=== IP-IN Pairs to Process ==="
for ITEM in "${IP_LIST[@]}"; do
    IFS=':' read -r IP_FILE REP <<< "$ITEM"
    BASENAME=$(basename "$IP_FILE" .bam)
    SAMPLE_NAME=$(echo "$BASENAME" | sed -E "s/_[Ii][Pp]${REP}.*//")
    
    shopt -s nullglob nocaseglob
    IN_CANDIDATES=("$BAM_DIR/${SAMPLE_NAME}_"[Ii][Nn]${REP}*.bam)
    shopt -u nullglob nocaseglob
    
    if [[ ${#IN_CANDIDATES[@]} -gt 0 ]]; then
        echo "  $SAMPLE_NAME rep$REP: $(basename "$IP_FILE") + $(basename "${IN_CANDIDATES[0]}")"
    else
        echo "  $SAMPLE_NAME rep$REP: $(basename "$IP_FILE") + [MISSING IN]" >&2
    fi
done
echo ""

# --------------------------
# Run MACS3 in parallel
# --------------------------
echo "Running MACS3 peak calling in parallel using $THREADS threads..."

if ! parallel -j "$THREADS" --halt soon,fail=1 --colsep ':' \
    run_macs3 {1} {2} ::: "${IP_LIST[@]}"; then
    echo "ERROR: Some MACS3 jobs failed. Check logs in $OUTPUT_DIR/logs" >&2
    exit 1
fi

# --------------------------
# Summary statistics
# --------------------------
echo ""
echo "=== MACS3 Peak Calling Summary ==="

shopt -s nullglob
PEAK_FILES=("$OUTPUT_DIR/macs3"/*_peaks.narrowPeak)
shopt -u nullglob

echo "Peak files generated: ${#PEAK_FILES[@]}"
echo ""

if [[ ${#PEAK_FILES[@]} -gt 0 ]]; then
    echo "=== Peak Counts per Sample ==="
    for peak_file in "${PEAK_FILES[@]}"; do
        sample=$(basename "$peak_file" _peaks.narrowPeak)
        count=$(wc -l < "$peak_file")
        echo "  $sample: $count peaks"
    done
fi

echo ""
echo "All MACS3 peak calling completed successfully."
echo "Results saved in: $OUTPUT_DIR/macs3"
echo "Logs saved in: $OUTPUT_DIR/logs"