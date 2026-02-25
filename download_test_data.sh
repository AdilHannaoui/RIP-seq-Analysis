#!/usr/bin/env bash
set -euo pipefail

# ============================================
# RNA-seq Test Dataset Downloader
# Source: ENA (EBI)
# ============================================

DATA_DIR="data/test"
BASE_URL="https://www.ebi.ac.uk/ena/browser/view"

mkdir -p "${DATA_DIR}"

download_sample () {
    local sample_name=$1
    local accession=$2
    
    echo "Downloading ${sample_name} (${accession})..."
    wget -q --show-progress \
        -O "${DATA_DIR}/${sample_name}.fastq.gz" \
        "${BASE_URL}/${accession}?download=true"
}

# Wild Type samples
declare -A WT=(
    [WT_1]=SRR3288156
    [WT_2]=SRR3288157
    [WT_3]=SRR3288167
    [WT_4]=SRR3288170
    [WT_5]=SRR3288171
    [WT_6]=SRR3288172
)

# Knockout samples
declare -A KO=(
    [KO_1]=SRR3288165
    [KO_2]=SRR3288159
    [KO_3]=SRR3288160
    [KO_4]=SRR3288161
)

for sample in "${!WT[@]}"; do
    download_sample "$sample" "${WT[$sample]}"
done

for sample in "${!KO[@]}"; do
    download_sample "$sample" "${KO[$sample]}"
done

echo "All downloads completed successfully."
