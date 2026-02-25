# ==========================
# RIP-seq Enrichment Analysis
# Author: Adil Hannaoui Anaaoui
# ==========================

library(dplyr)
library(readr)
library(purrr)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

message("=== [GO Enrichment] Starting ===")

# --------------------------
# Validate required variables
# --------------------------
if (!exists("WINDOWS_DIR") || is.null(WINDOWS_DIR)) {
  stop("ERROR: WINDOWS_DIR not defined in config.R")
}

if (!exists("DESEQ2_FC_DIR") || is.null(DESEQ2_FC_DIR)) {
  stop("ERROR: DESEQ2_FC_DIR not defined in config.R")
}

if (!exists("PLOTS_DIR") || is.null(PLOTS_DIR)) {
  stop("ERROR: PLOTS_DIR not defined in config.R")
}

if (!exists("GO_ONTOLOGY") || is.null(GO_ONTOLOGY)) {
  stop("ERROR: GO_ONTOLOGY not defined in config.R")
}

if (!exists("PVAL_CUTOFF") || is.null(PVAL_CUTOFF)) {
  stop("ERROR: PVAL_CUTOFF not defined in config.R")
}

if (!exists("QVAL_CUTOFF") || is.null(QVAL_CUTOFF)) {
  stop("ERROR: QVAL_CUTOFF not defined in config.R")
}

if (!exists("CONDITIONS") || is.null(CONDITIONS)) {
  stop("ERROR: CONDITIONS not defined in config.R")
}

# --------------------------
# Setup
# --------------------------
windows_dir <- WINDOWS_DIR
deseq_dir <- DESEQ2_FC_DIR

if (!dir.exists(windows_dir)) {
  stop("ERROR: WINDOWS_DIR not found: ", windows_dir)
}

if (!dir.exists(deseq_dir)) {
  stop("ERROR: DESEQ2_FC_DIR not found: ", deseq_dir)
}

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------
# Derive base conditions (LÓGICA ORIGINAL)
# --------------------------
ip_conditions <- CONDITIONS[grepl("_IP$", CONDITIONS)]
base_conditions <- unique(sub("_IP$", "", ip_conditions))

# --------------------------
# Function to run enrichment (LÓGICA ORIGINAL - SIN CAMBIOS)
# --------------------------
run_enrichment <- function(cond) {
  message("=== Enriquecimiento GO para condición: ", cond, " ===")
  
  # Archivo integrado filtrado (LÓGICA ORIGINAL)
  infile <- file.path(windows_dir, paste0("Integrated_", cond, "_filtered.tsv"))
  
  if (!file.exists(infile)) {
    warning("Integrated file not found: ", infile)
    return(NULL)
  }
  
  df <- tryCatch({
    read_tsv(infile, show_col_types = FALSE)
  }, error = function(e) {
    warning("Failed to read integrated file: ", e$message)
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  # Lista de genes (LÓGICA ORIGINAL)
  genes <- df$gene_id
  
  # Universo: todos los genes testeados por DESeq2 (LÓGICA ORIGINAL - EXACTAMENTE COMO ESTABA)
  deseq_file <- file.path(deseq_dir, paste0("DESeq2_", cond, "_IP_vs_IN_sig.csv"))
  
  if (!file.exists(deseq_file)) {
    warning("DESeq2 file not found: ", deseq_file)
    return(NULL)
  }
  
  deseq_all <- tryCatch({
    read_csv(deseq_file, show_col_types = FALSE)
  }, error = function(e) {
    warning("Failed to read DESeq2 file: ", e$message)
    return(NULL)
  })
  
  if (is.null(deseq_all)) return(NULL)
  
  universe <- deseq_all[[1]]   # primera columna = gene_id (LÓGICA ORIGINAL)
  
  # ==========================
  # Enriquecimiento GO (LÓGICA ORIGINAL - EXACTAMENTE COMO ESTABA)
  # ==========================
  
  ego <- tryCatch({
    enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = ORG_DB,
      keyType       = GENE_ID_TYPE,
      ont           = GO_ONTOLOGY,
      pAdjustMethod = "BH",
      pvalueCutoff  = PVAL_CUTOFF,
      qvalueCutoff  = QVAL_CUTOFF
    )
  }, error = function(e) {
    warning("GO enrichment failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(ego) || nrow(ego) == 0) {
    message("No hay términos GO para ", GO_ONTOLOGY, " en ", cond)
    return(NULL)
  }
  
  # Dotplot (LÓGICA ORIGINAL)
  p <- dotplot(ego, showCategory = 20) +
    ggtitle(paste0("GO ", GO_ONTOLOGY, " - ", cond))
  
  outfile <- file.path(
    PLOTS_DIR,
    paste0("GO_", GO_ONTOLOGY, "_", cond, "_dotplot.png")
  )
  
  ggsave(outfile, p, width = 8, height = 6)
  message("Dotplot generado: ", outfile)
  
  return(ego)
}

# --------------------------
# Run enrichment for all conditions (LÓGICA ORIGINAL)
# --------------------------
map(base_conditions, run_enrichment)


message("\n=== GO enrichment completed ===")

