#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Peak Intersection Module (parallel)
# Author: Adil Hannaoui Anaaoui
# ==========================

source "$(dirname "$0")/config.sh"

mkdir -p "$OUTPUT_DIR/macs2"
mkdir -p "$OUTPUT_DIR/logs"

cd "$WORKDIR"

#############################################
# Function to process a single sample
#############################################
process_sample() {
    BASENAME="$1"

    echo ">>> Processing sample: $BASENAME"

    # Collect narrowPeak files
    NARROWPEAK_FILES=()
    for NP in "$OUTPUT_DIR/macs2/${BASENAME}"_rep*_peaks.narrowPeak; do
        [[ -e "$NP" ]] || continue
        NARROWPEAK_FILES+=("$NP")
    done

    NUM_REPS=${#NARROWPEAK_FILES[@]}

    if [[ $NUM_REPS -lt 2 ]]; then
        echo "Not enough replicates for $BASENAME"
        return
    fi

    echo "Found $NUM_REPS replicates for $BASENAME"

    RAW="$OUTPUT_DIR/macs2/${BASENAME}_common_raw.bed"
    ANNOT="$OUTPUT_DIR/macs2/${BASENAME}_common_annotated.bed"
    COUNTS_TMP="$OUTPUT_DIR/macs2/${BASENAME}_common_counts.tmp"
    COUNTS_OUT="$OUTPUT_DIR/macs2/${BASENAME}_common_counts.txt"
    LOGFILE="$OUTPUT_DIR/logs/${BASENAME}_intersection.log"

    {
        #############################################
        # 3. CHAIN INTERSECT (rep1 ∩ rep2 ∩ rep3…)
        #############################################

        # Start with the first replicate
        cp "${NARROWPEAK_FILES[0]}" "$RAW"

        # Intersect sequentially with the rest
        for ((i=1; i<NUM_REPS; i++)); do
            bedtools intersect -a "$RAW" -b "${NARROWPEAK_FILES[$i]}" -wa -u > "${RAW}.tmp"
            mv "${RAW}.tmp" "$RAW"
        done

        # If RAW is empty → no common peaks
        if [[ ! -s "$RAW" ]]; then
            echo "No common peaks for $BASENAME"
            rm -f "$RAW"
            return
        fi

        #############################################
        # 4. Annotation
        #############################################
        bedtools intersect -a "$RAW" -b "$GTF_FILE" -wa -wb |
        awk -F'\t' '
        {
            gene_id="NA"; biotype="NA"
            if (match($0, /gene_id "([^"]+)"/, m)) gene_id=m[1]
            if (match($0, /gene_biotype "([^"]+)"/, b)) biotype=b[1]
            print $1, $2, $3, ".", ".", ".", gene_id, biotype
        }' OFS="\t" | sort -u > "$ANNOT"

        rm "$RAW"

        #############################################
        # 5. Header
        #############################################
        HEADER="chrom\tstart\tend\tname\tscore\tstrand\tGene_id\tBiotype"
        for REP in $(seq 1 "$NUM_REPS"); do
            HEADER+="\tIP${REP}\tIN${REP}"
        done

        #############################################
        # 6. Collect BAMs
        #############################################
        BAM_FILES=()
        for REP in $(seq 1 "$NUM_REPS"); do
            BAM_FILES+=("$OUTPUT_DIR/bowtie2/${BASENAME}_IP${REP}_trimmed.bam")
            BAM_FILES+=("$OUTPUT_DIR/bowtie2/${BASENAME}_IN${REP}_trimmed.bam")
        done

        #############################################
        # 7. multicov
        #############################################
        bedtools multicov -bams "${BAM_FILES[@]}" \
            -bed "$ANNOT" \
            > "$COUNTS_TMP"

        #############################################
        # 8. Combine header + counts
        #############################################
        {
            echo -e "$HEADER"
            cat "$COUNTS_TMP"
        } > "$COUNTS_OUT"

        rm "$COUNTS_TMP"

        echo ">>> Finished sample: $BASENAME"

    } > "$LOGFILE" 2>&1
}

export -f process_sample
export OUTPUT_DIR GTF_FILE

#############################################
# Detect all samples (rep1 only)
#############################################
SAMPLES=()

for FILE in "$OUTPUT_DIR/macs2"/*_rep1_peaks.narrowPeak; do
    [[ -e "$FILE" ]] || continue
    BASENAME=$(basename "$FILE" | sed 's/_rep1_peaks.narrowPeak//')
    SAMPLES+=("$BASENAME")
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "No samples found for intersection."
    exit 1
fi

echo "Found ${#SAMPLES[@]} samples."

#############################################
# Run in parallel
#############################################
parallel -j "$THREADS" process_sample ::: "${SAMPLES[@]}"

echo "All intersection analyses completed."
echo "Results saved in: $OUTPUT_DIR/macs2"
echo "Logs saved in: $OUTPUT_DIR/logs"
