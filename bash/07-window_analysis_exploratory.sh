#!/usr/bin/env bash
set -euo pipefail

# ==================================================
# RIP-seq Window Analysis (Exploratory)
# ==================================================

source "$(dirname "$0")/config.sh"

echo "=== [Window Analysis] Starting window-based RIP-seq enrichment ==="

WINDOW_DIR="${OUTPUT_DIR}/windows"
mkdir -p "${WINDOW_DIR}"

GENES_BED="${WINDOW_DIR}/genes.bed"
TSS_WINDOWS="${WINDOW_DIR}/genes_TSS_1kb.bed"

# ==================================================
# 1) Convert GTF → BED (gene-level)
# ==================================================
echo "[1/4] Extracting gene coordinates from GTF..."

awk '
    $3 == "gene" {
        # BED is 0-based start
        split($0, a, "\t")
        chr=a[1]
        start=a[4]-1
        end=a[5]
        gene_id="NA"

        # Extract gene_id from attributes
        match($0, /gene_id "([^"]+)"/, m)
        if (m[1] != "") gene_id=m[1]

        print chr"\t"start"\t"end"\t"gene_id
    }
' "$GTF_FILE" > "$GENES_BED"

echo "Gene BED written to: $GENES_BED"

# ==================================================
# 2) Create ±500 bp windows around TSS
# ==================================================
echo "[2/4] Creating ±500 bp TSS windows..."

# Generate genome sizes file if not present
GENOME_SIZES="${WINDOW_DIR}/genome.sizes"
cut -f1,2 "$BOWTIE2_INDEX.fa.fai" > "$GENOME_SIZES" 2>/dev/null || true

# If fai doesn't exist, warn
if [[ ! -f "$GENOME_SIZES" ]]; then
    echo "ERROR: genome.sizes not found. Ensure genome FASTA is indexed (.fai)."
    exit 1
fi

# Create windows
bedtools slop \
    -b 500 \
    -g "$GENOME_SIZES" \
    -i "$GENES_BED" \
    > "$TSS_WINDOWS"

echo "TSS windows written to: $TSS_WINDOWS"

# ==================================================
# 3) Compute coverage in windows for each merged IP BAM
# ==================================================
echo "[3/4] Computing coverage in windows..."

for bam in "$OUTPUT_DIR"/*_IP_merged_sorted.bam; do
    [[ -f "$bam" ]] || continue

    base=$(basename "$bam" .bam)
    out_cov="${WINDOW_DIR}/${base}_TSS_1kb_coverage.bed"

    echo "  - Processing $base"

    bedtools coverage \
        -a "$TSS_WINDOWS" \
        -b "$bam" \
        > "$out_cov"
done

echo "Coverage files written to: $WINDOW_DIR"

# ==================================================
# 4) Summary
# ==================================================
echo "[4/4] Window analysis completed."
echo "You can now integrate:"
echo "  - DESeq2 enriched genes"
echo "  - MACS2 peak-associated genes"
echo "  - Window-enriched genes"
echo "to identify high-confidence RIP-seq targets."
