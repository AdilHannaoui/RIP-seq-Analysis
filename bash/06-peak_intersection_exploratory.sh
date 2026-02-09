#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Peak Intersection Module (parallel)
# Author: Adil Hannaoui Anaaoui
# ==========================

source "$(dirname "$0")/config.sh"

# --------------------------
# Setup directories
# --------------------------
mkdir -p "$OUTPUT_DIR/macs3" "$OUTPUT_DIR/logs"
cd "$WORKDIR"

# --------------------------
# Validate dependencies
# --------------------------
for cmd in bedtools parallel; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "ERROR: Required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

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

MACS3_DIR="$OUTPUT_DIR/macs3"
BAM_DIR="$OUTPUT_DIR/bowtie2"

if [[ ! -d "$MACS3_DIR" ]]; then
    echo "ERROR: MACS3 directory not found at: $MACS3_DIR" >&2
    exit 1
fi

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found at: $BAM_DIR" >&2
    exit 1
fi

echo "Using GTF: $GTF"
echo "Using MACS3 directory: $MACS3_DIR"
echo "Using BAM directory: $BAM_DIR"

#############################################
# Function to process a single sample
#############################################
process_sample() {
    local BASENAME="$1"
    
    echo ">>> Processing sample: $BASENAME"
    
    # Collect narrowPeak files
    shopt -s nullglob
    local NARROWPEAK_FILES=("$MACS3_DIR/${BASENAME}"_rep*_peaks.narrowPeak)
    shopt -u nullglob
    
    local NUM_REPS=${#NARROWPEAK_FILES[@]}
    
    if [[ $NUM_REPS -eq 0 ]]; then
        echo "WARNING: No replicates found for $BASENAME" >&2
        return 0
    fi
    
    if [[ $NUM_REPS -lt 2 ]]; then
        echo "WARNING: Only $NUM_REPS replicate found for $BASENAME (need ≥2 for intersection)" >&2
        return 0
    fi
    
    echo "  Found $NUM_REPS replicates for $BASENAME"
    
    # Verify all peak files are not empty
    local EMPTY_COUNT=0
    for np in "${NARROWPEAK_FILES[@]}"; do
        if [[ ! -s "$np" ]]; then
            echo "  WARNING: Empty peak file: $(basename "$np")" >&2
            ((EMPTY_COUNT++))
        fi
    done
    
    if [[ $EMPTY_COUNT -eq $NUM_REPS ]]; then
        echo "  ERROR: All peak files are empty for $BASENAME" >&2
        return 1
    fi
    
    local RAW="$MACS3_DIR/${BASENAME}_common_raw.bed"
    local ANNOT="$MACS3_DIR/${BASENAME}_common_annotated.bed"
    local COUNTS_TMP="$MACS3_DIR/${BASENAME}_common_counts.tmp"
    local COUNTS_OUT="$MACS3_DIR/${BASENAME}_common_counts.txt"
    local LOGFILE="$OUTPUT_DIR/logs/${BASENAME}_intersection.log"
    
    {
        #############################################
        # Chain intersect (rep1 ∩ rep2 ∩ rep3…)
        #############################################
        
        # Start with the first replicate
        cp "${NARROWPEAK_FILES[0]}" "$RAW"
        
        # Intersect sequentially with the rest
        for ((i=1; i<NUM_REPS; i++)); do
            if [[ -s "${NARROWPEAK_FILES[$i]}" ]]; then
                bedtools intersect -a "$RAW" -b "${NARROWPEAK_FILES[$i]}" -wa -u > "${RAW}.tmp"
                mv "${RAW}.tmp" "$RAW"
            else
                echo "  WARNING: Skipping empty replicate: $(basename "${NARROWPEAK_FILES[$i]}")" >&2
            fi
        done
        
        # Check if RAW is empty → no common peaks
        if [[ ! -s "$RAW" ]]; then
            echo "  WARNING: No common peaks found for $BASENAME" >&2
            rm -f "$RAW"
            return 0
        fi
        
        local COMMON_PEAKS
        COMMON_PEAKS=$(wc -l < "$RAW")
        echo "  Found $COMMON_PEAKS common peaks"
        
        #############################################
        # Annotation
        #############################################
        if bedtools intersect -a "$RAW" -b "$GTF" -wa -wb | \
            awk -F'\t' '
            {
                gene_id="NA"; biotype="NA"
                if (match($0, /gene_id "([^"]+)"/, m)) gene_id=m[1]
                if (match($0, /gene_biotype "([^"]+)"/, b)) biotype=b[1]
                print $1, $2, $3, ".", ".", ".", gene_id, biotype
            }' OFS="\t" | sort -u > "$ANNOT"; then
            
            if [[ ! -s "$ANNOT" ]]; then
                echo "  WARNING: No annotated peaks for $BASENAME (peaks outside GTF regions)" >&2
            else
                local ANNOT_PEAKS
                ANNOT_PEAKS=$(wc -l < "$ANNOT")
                echo "  Annotated $ANNOT_PEAKS peaks"
            fi
        else
            echo "  ERROR: Annotation failed for $BASENAME" >&2
            rm -f "$RAW" "$ANNOT"
            return 1
        fi
        
        rm -f "$RAW"
        
        # If no annotated peaks, skip counting
        if [[ ! -s "$ANNOT" ]]; then
            rm -f "$ANNOT"
            return 0
        fi
        
        #############################################
        # Build header
        #############################################
        local HEADER="chrom\tstart\tend\tname\tscore\tstrand\tGene_id\tBiotype"
        for REP in $(seq 1 "$NUM_REPS"); do
            HEADER+="\tIP${REP}\tIN${REP}"
        done
        
        #############################################
        # Collect BAMs (case-insensitive, flexible naming)
        #############################################
        local BAM_FILES=()
        local MISSING_BAMS=0
        
        for REP in $(seq 1 "$NUM_REPS"); do
            # Try multiple naming patterns
            shopt -s nullglob nocaseglob
            local IP_CANDIDATES=("$BAM_DIR/${BASENAME}_"[Ii][Pp]${REP}*.bam)
            local IN_CANDIDATES=("$BAM_DIR/${BASENAME}_"[Ii][Nn]${REP}*.bam)
            shopt -u nullglob nocaseglob
            
            if [[ ${#IP_CANDIDATES[@]} -gt 0 ]]; then
                BAM_FILES+=("${IP_CANDIDATES[0]}")
            else
                echo "  ERROR: IP BAM not found for rep$REP: ${BASENAME}_IP${REP}*.bam" >&2
                ((MISSING_BAMS++))
                BAM_FILES+=("")  # Placeholder
            fi
            
            if [[ ${#IN_CANDIDATES[@]} -gt 0 ]]; then
                BAM_FILES+=("${IN_CANDIDATES[0]}")
            else
                echo "  ERROR: IN BAM not found for rep$REP: ${BASENAME}_IN${REP}*.bam" >&2
                ((MISSING_BAMS++))
                BAM_FILES+=("")  # Placeholder
            fi
        done
        
        if [[ $MISSING_BAMS -gt 0 ]]; then
            echo "  ERROR: Missing $MISSING_BAMS BAM file(s) for $BASENAME" >&2
            rm -f "$ANNOT"
            return 1
        fi
        
        # Verify all BAM files exist and are not empty
        for bam in "${BAM_FILES[@]}"; do
            if [[ ! -s "$bam" ]]; then
                echo "  ERROR: BAM file missing or empty: $bam" >&2
                rm -f "$ANNOT"
                return 1
            fi
        done
        
        #############################################
        # Run multicov
        #############################################
        if bedtools multicov -bams "${BAM_FILES[@]}" -bed "$ANNOT" > "$COUNTS_TMP"; then
            
            if [[ ! -s "$COUNTS_TMP" ]]; then
                echo "  ERROR: multicov produced empty output for $BASENAME" >&2
                rm -f "$ANNOT" "$COUNTS_TMP"
                return 1
            fi
            
            # Combine header + counts
            {
                echo -e "$HEADER"
                cat "$COUNTS_TMP"
            } > "$COUNTS_OUT"
            
            rm -f "$COUNTS_TMP" "$ANNOT"
            
            local FINAL_PEAKS
            FINAL_PEAKS=$(tail -n +2 "$COUNTS_OUT" | wc -l)
            echo "  Generated count matrix: $FINAL_PEAKS peaks"
            echo ">>> Finished sample: $BASENAME (OK)"
            
        else
            echo "  ERROR: multicov failed for $BASENAME" >&2
            rm -f "$ANNOT" "$COUNTS_TMP"
            return 1
        fi
        
    } > "$LOGFILE" 2>&1
    
    # Return success/failure based on log
    if grep -q "ERROR:" "$LOGFILE"; then
        return 1
    fi
    return 0
}

export -f process_sample
export OUTPUT_DIR GTF MACS3_DIR BAM_DIR

#############################################
# Detect all samples (rep1 as reference)
#############################################
echo "Searching for samples with replicates..."

shopt -s nullglob
REP1_FILES=("$MACS3_DIR"/*_rep1_peaks.narrowPeak)
shopt -u nullglob

SAMPLES=()
for FILE in "${REP1_FILES[@]}"; do
    BASENAME=$(basename "$FILE" _rep1_peaks.narrowPeak)
    SAMPLES+=("$BASENAME")
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: No samples with rep1 found for intersection in $MACS3_DIR" >&2
    echo "Searched for: *_rep1_peaks.narrowPeak" >&2
    exit 1
fi

echo "Found ${#SAMPLES[@]} samples with replicates."
echo ""

# Print sample summary
echo "=== Samples to Process ==="
for SAMPLE in "${SAMPLES[@]}"; do
    shopt -s nullglob
    REP_FILES=("$MACS3_DIR/${SAMPLE}"_rep*_peaks.narrowPeak)
    shopt -u nullglob
    echo "  $SAMPLE: ${#REP_FILES[@]} replicates"
done
echo ""

#############################################
# Run in parallel
#############################################
echo "Running peak intersection in parallel using $THREADS threads..."

if ! parallel -j "$THREADS" --halt soon,fail=1 process_sample ::: "${SAMPLES[@]}"; then
    echo "ERROR: Some intersection jobs failed. Check logs in $OUTPUT_DIR/logs" >&2
    exit 1
fi

# --------------------------
# Summary statistics
# --------------------------
echo ""
echo "=== Peak Intersection Summary ==="

shopt -s nullglob
COUNT_FILES=("$MACS3_DIR"/*_common_counts.txt)
shopt -u nullglob

echo "Common peak count files generated: ${#COUNT_FILES[@]}"

if [[ ${#COUNT_FILES[@]} -gt 0 ]]; then
    echo ""
    echo "=== Common Peaks per Sample ==="
    for count_file in "${COUNT_FILES[@]}"; do
        sample=$(basename "$count_file" _common_counts.txt)
        # Subtract header line
        peaks=$(($(wc -l < "$count_file") - 1))
        echo "  $sample: $peaks common peaks"
    done
fi

echo ""
echo "All intersection analyses completed successfully."
echo "Results saved in: $MACS3_DIR"
echo "Logs saved in: $OUTPUT_DIR/logs"