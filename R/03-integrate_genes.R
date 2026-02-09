# ==========================
# RIP-seq Data Integration
# Author: Adil Hannaoui Anaaoui
# ==========================

library(dplyr)
library(readr)
library(purrr)
library(tidyr)

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

message("=== [Data Integration] Starting ===")

# --------------------------
# Validate required variables
# --------------------------
if (!exists("CONDITIONS") || is.null(CONDITIONS) || length(CONDITIONS) == 0) {
  stop("ERROR: CONDITIONS not defined in config.R")
}

if (!exists("WINDOWS_DIR") || is.null(WINDOWS_DIR)) {
  stop("ERROR: WINDOWS_DIR not defined in config.R")
}

if (!exists("DESEQ2_FC_DIR") || is.null(DESEQ2_FC_DIR)) {
  stop("ERROR: DESEQ2_FC_DIR not defined in config.R")
}

if (!exists("DESEQ2_MACS3_DIR") || is.null(DESEQ2_MACS3_DIR)) {
  stop("ERROR: DESEQ2_MACS3_DIR not defined in config.R")
}

# --------------------------
# Setup
# --------------------------
windows_dir <- WINDOWS_DIR
deseq_dir   <- DESEQ2_FC_DIR
macs_dir    <- DESEQ2_MACS3_DIR

for (dir_path in list(windows_dir, deseq_dir, macs_dir)) {
  if (!dir.exists(dir_path)) {
    stop("ERROR: Directory not found: ", dir_path)
  }
}

dir.create(WINDOWS_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------
# Derive base conditions (LÓGICA ORIGINAL)
# --------------------------
ip_conditions <- CONDITIONS[grepl("_IP$", CONDITIONS)]
base_conditions <- unique(sub("_IP$", "", ip_conditions))
message("Conditions detected: ", paste(base_conditions, collapse = ", "))

# --------------------------
# Helper function (LÓGICA ORIGINAL)
# --------------------------
ensure_gene_id <- function(df) {
  first_col <- colnames(df)[1]
  df %>% rename(gene_id = !!first_col)
}

# --------------------------
# Function to process one condition (LÓGICA ORIGINAL)
# --------------------------
process_condition <- function(cond) {
  message("=== Procesando condición: ", cond, " ===")
  
  ## 1) DESeq2
  deseq_file <- file.path(deseq_dir, paste0("DESeq2_", cond, "_IP_vs_IN_sig.csv"))
  
  if (!file.exists(deseq_file)) {
    warning("DESeq2 file not found: ", deseq_file)
    return(NULL)
  }
  
  deseq <- tryCatch({
    read_csv(deseq_file) %>% ensure_gene_id()
  }, error = function(e) {
    warning("Failed to read DESeq2 file: ", e$message)
    return(NULL)
  })
  
  if (is.null(deseq)) return(NULL)
  
  deseq <- deseq %>%
    select(gene_id, log2FoldChange) %>%
    rename(
      !!paste0(cond, "_DESeq2_log2FC") := log2FoldChange
    )
  
  ## 2) MACS3
  macs_file <- file.path(macs_dir, paste0("DESeq2_res_", cond, "_sig.csv"))
  
  if (!file.exists(macs_file)) {
    warning("MACS3 file not found: ", macs_file)
    return(NULL)
  }
  
  macs <- tryCatch({
    read_csv(macs_file) %>% ensure_gene_id()
  }, error = function(e) {
    warning("Failed to read MACS3 file: ", e$message)
    return(NULL)
  })
  
  if (is.null(macs)) return(NULL)
  
  macs <- macs %>%
    select(gene_id, log2FoldChange) %>%
    rename(
      !!paste0(cond, "_MACS3_log2FC") := log2FoldChange
    )
  
  ## 3) TSS windows
  tss_file <- file.path(windows_dir, paste0(cond, "_IP_merged_sorted_TSS_1kb_coverage.bed"))
  
  if (!file.exists(tss_file)) {
    warning("TSS file not found: ", tss_file)
    return(NULL)
  }
  
  tss <- tryCatch({
    read_tsv(tss_file, col_names = FALSE)
  }, error = function(e) {
    warning("Failed to read TSS file: ", e$message)
    return(NULL)
  })
  
  if (is.null(tss)) return(NULL)
  
  tss <- tss %>%
    select(X4, X7) %>%
    rename(
      gene_id = X4,
      !!paste0(cond, "_TSS_coverage") := X7
    )
  
  ## 4) Integración (LÓGICA ORIGINAL)
  final <- deseq %>%
    full_join(macs, by = "gene_id") %>%
    full_join(tss,  by = "gene_id")
  
  ## 5) Eliminar genes incompletos (LÓGICA ORIGINAL)
  final <- final %>% drop_na()
  
  ## 6) Filtrar por log2FC de DESeq2 (LÓGICA ORIGINAL)
  logfc_col <- paste0(cond, "_DESeq2_log2FC")
  
  final <- final %>%
    filter(
      !!sym(logfc_col) > 1 |
        !!sym(logfc_col) < -1
    )
  
  ## 7) Ordenar por log2FC de DESeq2 (descendente) (LÓGICA ORIGINAL)
  final <- final %>%
    arrange(desc(!!sym(logfc_col)))
  
  ## 8) Guardar archivo final por condición (LÓGICA ORIGINAL)
  out_file <- file.path(windows_dir, paste0("Integrated_", cond, "_filtered.tsv"))
  write_tsv(final, out_file)
  message("Archivo final generado: ", out_file)
  
  return(final)
}

# --------------------------
# Process all conditions (LÓGICA ORIGINAL)
# --------------------------
all_results <- map(base_conditions, process_condition)

message("\n=== Data integration completed ===")