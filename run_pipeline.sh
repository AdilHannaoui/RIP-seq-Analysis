#!/usr/bin/env bash

set -e
set -o pipefail

mkdir -p logs

SECONDS=0
echo "=== Pipeline started successfully ===)"

echo "[1/11] Activating Conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ripseq-rpb4

echo "[2/11] Loading configs..."
CONFIG_R="R/config.R"
CONFIG_bash="bash/config.sh"

echo "[3/11] Running FastQC..."
bash bash/01-fastqc.sh $CONFIG_bash > logs/fastqc.log 2>&1

echo "[4/11] Running Trimmomatic..."
bash bash/02-trimming.sh $CONFIG_bash > logs/trimming.log 2>&1

echo "[5/11] Running Bowtie2..."
bash bash/03-alignment_bowtie2.sh $CONFIG_bash > logs/alignment.log 2>&1

echo "[6/11] Running FeatureCounts..."
bash bash/04-featurecounts.sh $CONFIG_bash > logs/featurecounts.log 2>&1

echo "[7/11] Running MACS3..."
bash bash/05-peak_calling_exploratory.sh $CONFIG_bash > logs/macs3.log 2>&1

echo "[8/11] Running Peak Intersection..."
bash bash/06-peak_intersection_exploratory.sh $CONFIG_bash > logs/intersection.log 2>&1

echo "[9/11] Running Window Analysis..."
bash bash/07-window_analysis_exploratory.sh $CONFIG_bash > logs/windows.log 2>&1

echo "[10/11] Running TSS metagene Analysis..."
bash bash/08-tss_metagene.sh $CONFIG_bash > logs/tss_metagene.log 2>&1

echo "[11/11] Running R analysis..."
Rscript R/01-load_counts.R $CONFIG_R > logs/load_counts.log 2>&1
Rscript R/02-deseq2_analysis.R $CONFIG_R > logs/deseq2.log 2>&1
Rscript R/02.1-deseq2_peaks_analysis.R $CONFIG_R > logs/deseq2_peaks.log 2>&1
Rscript R/03-integrate_genes.R $CONFIG_R > logs/integration.log 2>&1
Rscript R/04-enrichment_analysis.R $CONFIG_R > logs/enrichment.log 2>&1
Rscript R/05-visualizations.R $CONFIG_R > logs/visualizations.log 2>&1

duration=$SECONDS
echo "Total execution time: $((duration / 60)) minutes $((duration % 60)) seconds"

