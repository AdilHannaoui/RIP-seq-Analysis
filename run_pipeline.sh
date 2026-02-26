#!/usr/bin/env bash
set -euo pipefail

# ==========================
# RIP-seq Master Pipeline
# Author: Adil Hannaoui Anaaoui
# ==========================

PIPELINE_VERSION="v1.0.0"
START_TIME=$(date +%s)
LOG_DIR="logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ==========================
# Functions
# ==========================
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

run_step() {
    local step_num=$1
    local step_name=$2
    local script=$3
    local logfile="${LOG_DIR}/${step_num}-${step_name}_${TIMESTAMP}.log"
    
    log_info "[${step_num}/8] Running ${step_name}..."
    
    if bash "$script" > "$logfile" 2>&1; then
        log_info "✓ ${step_name} completed successfully"
        return 0
    else
        log_error "✗ ${step_name} failed. Check log: $logfile"
        return 1
    fi
}

run_r_step() {
    local step_num=$1
    local step_name=$2
    local script=$3
    local logfile="${LOG_DIR}/${step_num}-${step_name}_${TIMESTAMP}.log"
    
    log_info "[${step_num}/16] Running ${step_name}..."
    
    if Rscript "$script" > "$logfile" 2>&1; then
        log_info "✓ ${step_name} completed successfully"
        return 0
    else
        log_error "✗ ${step_name} failed. Check log: $logfile"
        return 1
    fi
}

# ==========================
# Main Pipeline
# ==========================
echo "=========================================="
echo "  RIP-seq Analysis Pipeline ${PIPELINE_VERSION}"
echo "  Started: $(date)"
echo "=========================================="

# Create logs directory
mkdir -p "$LOG_DIR"

# 1. Environment Activation
log_info "[1/16] Activating Conda environment..."
if ! source ~/miniconda3/etc/profile.d/conda.sh; then
    log_error "Failed to source conda.sh"
    exit 1
fi

if ! conda activate ripseq-rpb4; then
    log_error "Failed to activate conda environment"
    exit 1
fi

log_info "Environment activated: $(conda info --envs | grep '*' | awk '{print $1}')"

# 2. Validate configuration files
log_info "[2/16] Validating configuration files..."
if [[ ! -f "bash/config.sh" ]]; then
    log_error "bash/config.sh not found"
    exit 1
fi

if [[ ! -f "R/config.R" ]]; then
    log_error "R/config.R not found"
    exit 1
fi

log_info "Configuration files validated"

# 3-6. Bash steps
run_step "3" "fastqc" "bash/01-fastqc.sh" || exit 1
run_step "4" "trimming" "bash/02-trimming.sh" || exit 1
run_step "5" "alignment" "bash/03-alignment_bowtie2.sh" || exit 1
run_step "6" "featurecounts" "bash/04-featurecounts.sh" || exit 1
run_step "7" "featurecounts" "bash/05-peak_calling_exploratory.sh" || exit 1
run_step "8" "featurecounts" "bash/06-peak_intersection_exploratory.sh" || exit 1
run_step "9" "featurecounts" "bash/07-window_analysis_exploratory.sh" || exit 1
run_step "10" "featurecounts" "bash/08-tss_metagene.sh" || exit 1

# 11-16. R steps
run_r_step "11" "load_counts" "R/01-load_counts.R" || exit 1
run_r_step "12" "deseq2" "R/02-deseq2_analysis.R" || exit 1
run_r_step "13" "deseq2" "R/02.1-deseq2_peaks_analysis.R" || exit 1
run_r_step "14" "enrichment" "R/03-integrate_genes.R" || exit 1
run_r_step "15" "enrichment" "R/04-enrichment_analysis.R" || exit 1
run_r_step "16" "visualization" "R/05-visualizations.R" || exit 1

# ==========================
# Completion
# ==========================
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "  ✓ Pipeline completed successfully!"
echo "  Total time: $((DURATION / 60))m $((DURATION % 60))s"
echo "  Finished: $(date)"
echo "=========================================="
echo ""
log_info "Logs saved in: $LOG_DIR/"
log_info "Results available in: Plots/ and data/" seconds"
