# ==========================
# RIP-seq Analysis Config
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Project structure
# --------------------------
PROJECT_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))

DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")

COUNTS_MATRIX_PATH <- file.path(OUTPUT_DIR, "counts_matrix.rds")
SAMPLE_METADATA_PATH <- file.path(OUTPUT_DIR, "colData.rds")
FEATURECOUNTS_DIR <- file.path(OUTPUT_DIR, "featurecounts")
PLOTS_DIR <- file.path(PROJECT_ROOT, "Plots")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# --------------------------
# Experimental design
# --------------------------
CONDITIONS <- c(
  rep("Rpb4-S/T-A", 3),
  rep("Rpb4-S/T-D", 3),
  rep("WT", 6)
)

REFERENCE_CONDITION <- "WT"

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
