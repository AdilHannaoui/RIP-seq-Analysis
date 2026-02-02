#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq MACS2 Peak Calling (parallel)
# Author: Adil Hannaoui Anaaoui
# ==========================

# Load global config
source "$(dirname "$0")/config.sh"

mkdir -p "$OUTPUT_DIR/macs2"
mkdir -p "$OUTPUT_DIR/logs"

cd "$WORKDIR"

# --------------------------
# Function to run MACS2 for a single IP/IN pair
# --------------------------
run_macs2() {
    IP_FILE="$1"
    REP="$2"

    BASENAME=$(basename "$IP_FILE")
    SAMPLE_NAME=$(echo "$BASENAME" | sed -E "s/_[Ii][PNpn]${REP}.*//")

    # Find matching IN file
    IN_FILE=$(ls "$OUTPUT_DIR/bowtie2"/"${SAMPLE_NAME}"_[Ii][Nn]${REP}*.bam 2>/dev/null || true)

    if [[ ! -f "$IN_FILE" ]]; then
        echo "ERROR: No IN file found for sample $SAMPLE_NAME (rep $REP)"
        return
    fi

    LOGFILE="$OUTPUT_DIR/logs/macs2_${SAMPLE_NAME}_rep${REP}.log"

    echo ">>> Running MACS2 for $SAMPLE_NAME (rep $REP)"

    macs2 callpeak \
        -t "$IP_FILE" \
        -c "$IN_FILE" \
        --format BAM \
        --name "$OUTPUT_DIR/macs2/${SAMPLE_NAME}_rep${REP}" \
        --pvalue 0.1 \
        --call-summits \
        --nomodel \
        --extsize 150 \
        --gsize 1.2e7 \
        > "$LOGFILE" 2>&1

    echo ">>> Finished MACS2 for $SAMPLE_NAME (rep $REP)"
}

export -f run_macs2
export OUTPUT_DIR

# --------------------------
# Build list of IP files for all replicates
# --------------------------
IP_LIST=()

for REP in $(seq 1 "$MAX_REPLICAS"); do
    for IP_FILE in "$OUTPUT_DIR/bowtie2"/*_[Ii][Pp]${REP}*.bam; do
        [[ -e "$IP_FILE" ]] || continue
        IP_LIST+=("$IP_FILE:$REP")
    done
done

if [[ ${#IP_LIST[@]} -eq 0 ]]; then
    echo "No IP files found for MACS2 peak calling."
    exit 1
fi

echo "Found ${#IP_LIST[@]} IP files for MACS2."

# --------------------------
# Run MACS2 in parallel
# --------------------------
echo "Running MACS2 peak calling in parallel using $THREADS threads..."

parallel -j "$THREADS" --colsep ':' \
    run_macs2 {1} {2} ::: "${IP_LIST[@]}"


echo "All MACS2 peak calling completed."
echo "Results saved in: $OUTPUT_DIR/macs2"
echo "Logs saved in: $OUTPUT_DIR/logs"
