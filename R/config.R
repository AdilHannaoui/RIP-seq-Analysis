# ==========================
# RIP-seq Analysis Config
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Project structure
# --------------------------
PROJECT_ROOT <- getwd()

DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
COUNTS_FILES <- list.files(
    path = "output/macs2",
    pattern = "_common_counts\\.txt$",
    full.names = TRUE
)
SAMPLE_METADATA_PATH <- "output/colData.rds"
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")

CONDITIONS <- c(
  rep("M1_IN", 3),
  rep("M1_IP", 3),
  rep("M12_IN", 3),
  rep("M12_IP", 3),
  rep("WT_IN", 3),
  rep("WT_IP", 3)
)

REFERENCE_CONDITION <- "IN"

library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

REFERENCE_GENOME <- BSgenome.Scerevisiae.UCSC.sacCer3
REFERENCE_TXDB <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

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
PADJ_METHOD <- "BH"

# --------------------------
# Organism database
# --------------------------
organism <- "yeast"
ORG_DB <- "org.Sc.sgd.db"
GENE_ID_TYPE <- "ORF" 
KEGG_ORGANISM <- "sce"
