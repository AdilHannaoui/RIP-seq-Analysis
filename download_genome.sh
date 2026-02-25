#!/usr/bin/env bash
set -euo pipefail
# ============================================
# Reference Genome Downloader
# Source: Ensembl release 112
# Organism: Saccharomyces cerevisiae R64-1-1
# ============================================

REF_DIR="Bowtie2/cerevisiae"
BASE_URL="https://ftp.ensembl.org/pub/release-112"
FASTA_GZ="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
GTF_GZ="Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz"

mkdir -p "${REF_DIR}"

download_reference() {
    local filename=$1
    local url=$2

    echo "Downloading ${filename}..."
    wget -q --show-progress \
        -O "${REF_DIR}/${filename}" \
        "${url}"
}

# Download FASTA
download_reference "${FASTA_GZ}" \
    "${BASE_URL}/fasta/saccharomyces_cerevisiae/dna/${FASTA_GZ}"

# Download GTF
download_reference "${GTF_GZ}" \
    "${BASE_URL}/gtf/saccharomyces_cerevisiae/${GTF_GZ}"

# Decompress
echo "Decompressing files..."
gunzip "${REF_DIR}/${FASTA_GZ}"
gunzip "${REF_DIR}/${GTF_GZ}"

echo "Building Bowtie2 index..."
bowtie2-build \
    "${REF_DIR}/${FASTA_GZ%.gz}" \
    "${REF_DIR}/index/genome"

echo "Reference genome ready in: ${REF_DIR}"
