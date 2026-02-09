# ==========================
# RIP-seq Analysis Config
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Project structure
# --------------------------
PROJECT_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))

DATA_DIR <- file.path(PROJECT_ROOT, "data")
GTF_DIR <- file.path(PROJECT_ROOT, "Bowtie2")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
DESEQ2_FC_DIR <- file.path(OUTPUT_DIR, "DESeq2 FeatureCounts")
BAM_DIR <- file.path(OUTPUT_DIR, "bowtie2")
DESEQ2_MACS3_DIR <- file.path(OUTPUT_DIR, "DESeq2 Macs3")
WINDOWS_DIR <- file.path(OUTPUT_DIR, "Windows Analysis")


GTF <- file.path(GTF_DIR, "cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf")
COUNTS_MATRIX_PATH_FC <- file.path(DESEQ2_FC_DIR, "counts_matrix.rds")
SAMPLE_METADATA_PATH_FC <- file.path(DESEQ2_FC_DIR, "colData.rds")
DDS_FC <- file.path(DESEQ2_FC_DIR, "dds.rds")
FEATURECOUNTS_DIR <- file.path(OUTPUT_DIR, "featurecounts")
PLOTS_DIR <- file.path(OUTPUT_DIR, "Plots")

# --------------------------
# Experimental design
# --------------------------
CONDITIONS <- c(
  rep("M1_IN", 3),
  rep("M1_IP", 3),
  rep("M12_IN", 3),
  rep("M12_IP", 3),
  rep("WT_IN", 3),
  rep("WT_IP", 3)
)

REFERENCE_CONDITION <- "IN"

# --------------------------
# DESeq2 parameters
# --------------------------
PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1
MIN_COUNTS_FILTER <- 10

# --------------------------
# Enrichment analysis
# --------------------------
GO_ONTOLOGY <- "BP" 
PVAL_CUTOFF <- 0.05
QVAL_CUTOFF <- 0.05

# --------------------------
# Organism database
# --------------------------
ORG_DB <- "org.Sc.sgd.db"
GENE_ID_TYPE <- "ORF" 
