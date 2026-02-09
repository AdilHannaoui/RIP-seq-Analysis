#!/bin/bash

# ==================================================
# Global configuration file for RIP-seq pipeline
# ==================================================

# ----------------------
# Project directories
# ----------------------
WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTQ_DIR="$WORKDIR/data"
OUTPUT_DIR="$WORKDIR/output"
BAM_DIR="$OUTPUT_DIR/bowtie2"
PLOTS_DIR="$OUTPUT_DIR/Plots"

# ----------------------
# Computational resources
# ----------------------
THREADS=10
MAX_REPLICAS=3

# ----------------------
# Tools and references
# ----------------------

# Trimmomatic
TRIMMO_JAR="$CONDA_PREFIX/share/trimmomatic/trimmomatic.jar"
ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-SE.fa"

# BOWTIE2
BOWTIE2_INDEX="$WORKDIR/Bowtie2/cerevisiae/index/genome"

# Gene annotation
GTF="$WORKDIR/Bowtie2/cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf"
GENOME_FA="$WORKDIR/Bowtie2/cerevisiae/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"

# ----------------------
# Output subdirectories
# ----------------------
FASTQC_PRE_DIR="$OUTPUT_DIR/fastqc_pre"
FASTQC_POST_DIR="$OUTPUT_DIR/fastqc_post"
FASTQ_TRIM="$OUTPUT_DIR/fastq_trimmed"
BOWTIE2_DIR="$OUTPUT_DIR/bowtie2"
LOG_DIR="$OUTPUT_DIR/logs"
