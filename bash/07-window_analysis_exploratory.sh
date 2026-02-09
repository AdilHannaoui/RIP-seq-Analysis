#!/usr/bin/env bash
set -euo pipefail

# ==================================================
# RIP-seq Window Analysis (TSS-centered, 1 kb windows)
# With automatic merge of IP replicates per condition
# ==================================================

source "$(dirname "$0")/config.sh"

echo "=== [Window Analysis] Starting TSS-centered RIP-seq enrichment ==="

# --------------------------
# Setup directories
# --------------------------
WINDOW_DIR="${OUTPUT_DIR}/Windows Analysis"
BAM_DIR="${OUTPUT_DIR}/bowtie2"
mkdir -p "$WINDOW_DIR"

# --------------------------
# Validate dependencies
# --------------------------
for cmd in samtools bedtools awk; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "ERROR: Required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# --------------------------
# Validate required files
# --------------------------
if [[ -z "${GTF:-}" ]]; then
    echo "ERROR: GTF not defined in config.sh" >&2
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF file not found at: $GTF" >&2
    exit 1
fi

if [[ -z "${GENOME_FA:-}" ]]; then
    echo "ERROR: GENOME_FA not defined in config.sh" >&2
    exit 1
fi

if [[ ! -f "$GENOME_FA" ]]; then
    echo "ERROR: Genome FASTA not found at: $GENOME_FA" >&2
    exit 1
fi

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory not found at: $BAM_DIR" >&2
    exit 1
fi

echo "Using GTF: $GTF"
echo "Using Genome: $GENOME_FA"
echo "Using BAM directory: $BAM_DIR"
echo ""

# Define output files
GENES_TSS="${WINDOW_DIR}/genes_TSS.bed"
TSS_WINDOWS="${WINDOW_DIR}/genes_TSS_1kb.bed"
GENOME_SIZES="${WINDOW_DIR}/genome.sizes"

# ==================================================
# 0) Auto-detect conditions and merge IP replicates
# ==================================================
echo "[0/4] Detecting conditions and merging replicates..."

cd "$BAM_DIR"

# Auto-detect conditions from BAM files (case-insensitive)
shopt -s nullglob nocaseglob
IP_BAMS=(*_[Ii][Pp]*.bam)
shopt -u nullglob nocaseglob

if [[ ${#IP_BAMS[@]} -eq 0 ]]; then
    echo "ERROR: No IP BAM files found in $BAM_DIR" >&2
    exit 1
fi

# Extract unique conditions
declare -A conditions_map
for bam in "${IP_BAMS[@]}"; do
    # Extract condition name (everything before _IP)
    cond=$(echo "$bam" | sed -E 's/_[Ii][Pp].*//')
    conditions_map["$cond"]=1
done

# Convert to array
conditions=("${!conditions_map[@]}")

if [[ ${#conditions[@]} -eq 0 ]]; then
    echo "ERROR: No conditions detected from BAM files" >&2
    exit 1
fi

echo "Detected ${#conditions[@]} conditions: ${conditions[*]}"
echo ""

# --------------------------
# Function to merge BAMs (IP or IN)
# --------------------------
merge_bams() {
    local cond="$1"
    local type="$2"  # IP or IN
    
    # Find replicates (case-insensitive)
    shopt -s nullglob nocaseglob
    local reps=("${cond}_${type}"*.bam)
    shopt -u nullglob nocaseglob
    
    local merged_bam="${cond}_${type}_merged_sorted.bam"
    
    if [[ ${#reps[@]} -eq 0 ]]; then
        echo "  WARNING: No ${type} BAMs found for $cond" >&2
        return 1
    elif [[ ${#reps[@]} -eq 1 ]]; then
        echo "  ${type}: 1 replicate → symlinking as merged"
        ln -sf "$(basename "${reps[0]}")" "$merged_bam"
        
        # Index if not exists
        if [[ ! -f "${merged_bam}.bai" ]]; then
            if ! samtools index "$merged_bam"; then
                echo "  ERROR: Failed to index $merged_bam" >&2
                return 1
            fi
        fi
    else
        echo "  ${type}: ${#reps[@]} replicates → merging"
        
        # Verify all BAMs are not empty
        for bam in "${reps[@]}"; do
            if [[ ! -s "$bam" ]]; then
                echo "  ERROR: Empty BAM file: $bam" >&2
                return 1
            fi
        done
        
        local temp_merged="${cond}_${type}_merged.bam"
        
        # Merge
        if ! samtools merge "$temp_merged" "${reps[@]}"; then
            echo "  ERROR: Failed to merge ${type} BAMs for $cond" >&2
            rm -f "$temp_merged"
            return 1
        fi
        
        # Sort
        if ! samtools sort -@ "${THREADS:-4}" -o "$merged_bam" "$temp_merged"; then
            echo "  ERROR: Failed to sort merged ${type} BAM for $cond" >&2
            rm -f "$temp_merged" "$merged_bam"
            return 1
        fi
        
        rm -f "$temp_merged"
        
        # Index
        if ! samtools index "$merged_bam"; then
            echo "  ERROR: Failed to index $merged_bam" >&2
            return 1
        fi
        
        echo "  Done: $merged_bam"
    fi
    
    return 0
}

export -f merge_bams
export THREADS

# Process each condition
for cond in "${conditions[@]}"; do
    echo "Processing condition: $cond"
    
    # Merge IP
    if ! merge_bams "$cond" "IP"; then
        echo "  WARNING: Failed to process IP for $cond" >&2
    fi
    
    # Merge IN
    if ! merge_bams "$cond" "IN"; then
        echo "  WARNING: Failed to process IN for $cond" >&2
    fi
    
    echo ""
done

cd "$WORKDIR"

# Verify at least one merged IP BAM was created
shopt -s nullglob
MERGED_IP_BAMS=("$BAM_DIR"/*_IP_merged_sorted.bam)
shopt -u nullglob

if [[ ${#MERGED_IP_BAMS[@]} -eq 0 ]]; then
    echo "ERROR: No merged IP BAMs were created" >&2
    exit 1
fi

echo "Successfully created ${#MERGED_IP_BAMS[@]} merged IP BAM(s)"
echo ""

# ==================================================
# 1) Extract TSS from GTF (strand-aware)
# ==================================================
echo "[1/4] Extracting TSS coordinates from GTF..."

if awk '
  $3=="gene" {
    chr=$1; start=$4; end=$5; strand=$7; gid="NA";
    match($0, /gene_id "([^"]+)"/, m); if (m[1]!="") gid=m[1];
    
    if (strand=="+") { tss = start - 1 }
    else { tss = end - 1 }
    
    # Filter out invalid chromosomes/coordinates
    if (tss >= 0 && chr !~ /^#/) {
        print chr"\t"tss"\t"tss+1"\t"gid"\t0\t"strand
    }
  }
' "$GTF" > "$GENES_TSS"; then
    
    if [[ ! -s "$GENES_TSS" ]]; then
        echo "ERROR: No TSS coordinates extracted from GTF" >&2
        exit 1
    fi
    
    TSS_COUNT=$(wc -l < "$GENES_TSS")
    echo "Extracted $TSS_COUNT TSS coordinates"
    echo "TSS BED written to: $GENES_TSS"
else
    echo "ERROR: Failed to extract TSS from GTF" >&2
    exit 1
fi

echo ""

# ==================================================
# 2) Generate genome sizes
# ==================================================
echo "[2/4] Generating genome sizes..."

if [[ ! -f "${GENOME_FA}.fai" ]]; then
    echo "FASTA index not found. Creating ${GENOME_FA}.fai ..."
    if ! samtools faidx "${GENOME_FA}"; then
        echo "ERROR: Failed to index genome FASTA" >&2
        exit 1
    fi
fi

if cut -f1,2 "${GENOME_FA}.fai" > "$GENOME_SIZES"; then
    if [[ ! -s "$GENOME_SIZES" ]]; then
        echo "ERROR: Empty genome sizes file" >&2
        exit 1
    fi
    echo "Genome sizes written to: $GENOME_SIZES"
else
    echo "ERROR: Failed to generate genome sizes" >&2
    exit 1
fi

echo ""

# ==================================================
# 3) Create ±500 bp windows around TSS
# ==================================================
echo "[3/4] Creating ±500 bp TSS windows..."

if bedtools slop \
    -b 500 \
    -g "$GENOME_SIZES" \
    -i "$GENES_TSS" \
    > "$TSS_WINDOWS"; then
    
    if [[ ! -s "$TSS_WINDOWS" ]]; then
        echo "ERROR: Empty TSS windows file" >&2
        exit 1
    fi
    
    WINDOW_COUNT=$(wc -l < "$TSS_WINDOWS")
    echo "Created $WINDOW_COUNT TSS windows (±500 bp)"
    echo "TSS windows written to: $TSS_WINDOWS"
else
    echo "ERROR: Failed to create TSS windows" >&2
    exit 1
fi

echo ""

# ==================================================
# 4) Compute coverage for MERGED IP BAMs
# ==================================================
echo "[4/4] Computing coverage in windows (merged IP only)..."

# Function to compute coverage for a single BAM
compute_coverage() {
    local bam="$1"
    local base
    base=$(basename "$bam" .bam)
    
    local out_cov="${WINDOW_DIR}/${base}_TSS_1kb_coverage.bed"
    
    echo "  Processing $base"
    
    # Verify BAM is indexed
    if [[ ! -f "${bam}.bai" ]]; then
        echo "    WARNING: Index missing, creating it..." >&2
        if ! samtools index "$bam"; then
            echo "    ERROR: Failed to index $bam" >&2
            return 1
        fi
    fi
    
    # Compute coverage
    if bedtools coverage \
        -a "$TSS_WINDOWS" \
        -b "$bam" \
        -counts \
        > "$out_cov"; then
        
        if [[ ! -s "$out_cov" ]]; then
            echo "    WARNING: Empty coverage file for $base" >&2
            return 1
        fi
        
        local total_counts
        total_counts=$(awk '{sum+=$NF} END {print sum}' "$out_cov")
        echo "    Done: $total_counts total counts"
        return 0
    else
        echo "    ERROR: Coverage computation failed for $base" >&2
        return 1
    fi
}

export -f compute_coverage
export WINDOW_DIR TSS_WINDOWS

# Run coverage computation
COVERAGE_FAILED=0

for bam in "${MERGED_IP_BAMS[@]}"; do
    if ! compute_coverage "$bam"; then
        ((COVERAGE_FAILED++))
    fi
done

echo ""

# --------------------------
# Summary
# --------------------------
echo "=== Window Analysis Summary ==="

shopt -s nullglob
COVERAGE_FILES=("$WINDOW_DIR"/*_TSS_1kb_coverage.bed)
shopt -u nullglob

echo "Coverage files generated: ${#COVERAGE_FILES[@]}"

if [[ $COVERAGE_FAILED -gt 0 ]]; then
    echo "WARNING: $COVERAGE_FAILED coverage computation(s) failed" >&2
fi

if [[ ${#COVERAGE_FILES[@]} -gt 0 ]]; then
    echo ""
    echo "=== Coverage Statistics ==="
    for cov_file in "${COVERAGE_FILES[@]}"; do
        sample=$(basename "$cov_file" _TSS_1kb_coverage.bed)
        regions=$(wc -l < "$cov_file")
        total_counts=$(awk '{sum+=$NF} END {print sum}' "$cov_file")
        avg_counts=$(awk -v tot="$total_counts" -v reg="$regions" 'BEGIN {printf "%.2f", tot/reg}')
        echo "  $sample: $regions regions, $total_counts total counts (avg: $avg_counts)"
    done
fi

echo ""
echo "All output files in: $WINDOW_DIR"
echo "=== Window analysis completed ==="

# Return error if all coverage computations failed
if [[ ${#COVERAGE_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No coverage files were generated" >&2
    exit 1
fi