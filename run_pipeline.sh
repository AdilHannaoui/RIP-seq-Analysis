#!/usr/bin/env bash

set -e
set -o pipefail

SECONDS=0
echo "=== Pipeline started successfully ===)"

# 1. Environment Activation
echo "[1/10] Activating Conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ripseq-rpb4

# 2. Configuration loading
echo "[2/10] Loading R config..."
CONFIG_R="R/config.R"
echo "[2/7] Loading bash config..."
CONFIG_bash="bash/config.sh"


# 3. QC
echo "[3/10] Running FastQC..."
bash bash/01-fastqc.sh $CONFIG_bash

# 4. Trimming
echo "[4/10] Running Trimmomatic..."
bash bash/02-trimming.sh $CONFIG_bash

# 5. Alignment
echo "[5/10] Running Bowtie2..."
bash bash/03-alignment_bowtie2.sh $CONFIG_bash

# 6. Peak Calling
echo "[6/10] Running macs2..."
bash bash/04-featurecounts.sh $CONFIG_bash

# 7. Peak Intersection
echo "[7/10] Making Peak Intersection..."
bash bash/05-peak_calling_exploratory.sh $CONFIG_bash

# 8. R Analysis
echo "[8/10] Running DESeq2 and downstream analysis..."
Rscript R/01-load_counts.R $CONFIG_R
Rscript R/02-deseq2_analysis.R $CONFIG_R
Rscript R/02.1-deseq2_peaks_analysis.R $CONFIG_R


duration=$SECONDS
echo "Total execution time: $((duration / 60)) minutes $((duration % 60)) seconds"
