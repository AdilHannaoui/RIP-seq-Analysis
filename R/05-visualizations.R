# ==========================
# RIP-seq Visualization: PCA, Heatmap, Upset
# Author: Adil Hannaoui Anaaoui
# ==========================

library(DESeq2)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(pheatmap)
library(UpSetR)

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

message("=== [Visualization] Starting ===")

# --------------------------
# Validate required variables
# --------------------------
if (!exists("PLOTS_DIR") || is.null(PLOTS_DIR)) {
  stop("ERROR: PLOTS_DIR not defined in config.R")
}

if (!exists("WINDOWS_DIR") || is.null(WINDOWS_DIR)) {
  stop("ERROR: WINDOWS_DIR not defined in config.R")
}

if (!exists("COUNTS_MATRIX_PATH_FC") || is.null(COUNTS_MATRIX_PATH_FC)) {
  stop("ERROR: COUNTS_MATRIX_PATH_FC not defined in config.R")
}

if (!exists("SAMPLE_METADATA_PATH_FC") || is.null(SAMPLE_METADATA_PATH_FC)) {
  stop("ERROR: SAMPLE_METADATA_PATH_FC not defined in config.R")
}

if (!exists("DDS_FC") || is.null(DDS_FC)) {
  stop("ERROR: DDS_FC not defined in config.R")
}

# Create output directory
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------
# Validate files exist
# --------------------------
if (!file.exists(COUNTS_MATRIX_PATH_FC)) {
  stop("ERROR: Counts matrix not found at: ", COUNTS_MATRIX_PATH_FC)
}

if (!file.exists(SAMPLE_METADATA_PATH_FC)) {
  stop("ERROR: Sample metadata not found at: ", SAMPLE_METADATA_PATH_FC)
}

if (!file.exists(DDS_FC)) {
  stop("ERROR: DESeq2 object not found at: ", DDS_FC)
}

if (!dir.exists(WINDOWS_DIR)) {
  stop("ERROR: Windows directory not found at: ", WINDOWS_DIR)
}

# ==========================
# Helper function: generate PCA
# ==========================
generate_pca <- function(vsd, colData, subset_pattern = NULL, filename, width = 7, height = 6) {
  
  if (!is.null(subset_pattern)) {
    message(sprintf("Generating PCA (%s only)...", subset_pattern))
    vsd_subset <- vsd[, grepl(subset_pattern, colData$condition)]
  } else {
    message("Generating global PCA...")
    vsd_subset <- vsd
  }
  
  pca_plot <- plotPCA(vsd_subset, intgroup = "condition")
  outfile <- file.path(PLOTS_DIR, filename)
  
  tryCatch({
    ggsave(outfile, pca_plot, width = width, height = height)
    message("PCA saved: ", outfile)
  }, error = function(e) {
    warning("Failed to save PCA plot: ", e$message)
  })
  
  return(pca_plot)
}

# ==========================
# Helper function: read integrated TSVs
# ==========================
read_integrated_tsvs <- function(dir_path, pattern = "Integrated_.*_filtered\\.tsv$") {
  tsv_files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  
  if (length(tsv_files) == 0) {
    stop("No integrated TSV files found in ", dir_path)
  }
  
  return(tsv_files)
}

# ==========================
# Load data
# ==========================
message("\n[1/4] Loading data...")

countData <- readRDS(COUNTS_MATRIX_PATH_FC)
colData   <- readRDS(SAMPLE_METADATA_PATH_FC)
dds       <- readRDS(DDS_FC)

# Validate data
if (nrow(countData) == 0 || ncol(countData) == 0) {
  stop("ERROR: Empty counts matrix")
}

if (nrow(colData) == 0) {
  stop("ERROR: Empty sample metadata")
}

message("Loaded ", nrow(countData), " genes x ", ncol(countData), " samples")

# Setup sample metadata (LÓGICA ORIGINAL)
colData$sample <- colnames(countData)
rownames(colData) <- colData$sample
colData$condition <- gsub("[0-9]+$", "", colData$sample)

# VST transformation
message("Performing variance stabilizing transformation...")
vsd <- vst(dds, blind = FALSE)

# ==========================
# Generate PCA plots
# ==========================
message("\n[2/4] === PCA Analysis ===")

generate_pca(vsd, colData, "IP", "PCA_IP_only.png", width = 7, height = 6)
generate_pca(vsd, colData, "IN", "PCA_IN_only.jpg", width = 10, height = 8)

# ==========================
# Integrated Heatmap
# ==========================
message("\n[3/4] === Integrated Heatmap ===")

tsv_files <- read_integrated_tsvs(WINDOWS_DIR)
message("Found ", length(tsv_files), " integrated TSV files")

# Read and process TSVs
df_list <- lapply(tsv_files, function(f) {
  cond <- sub("Integrated_(.*)_filtered\\.tsv", "\\1", basename(f))
  
  df <- tryCatch({
    read_tsv(f, show_col_types = FALSE)
  }, error = function(e) {
    warning("Failed to read ", basename(f), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(df) || nrow(df) == 0) {
    warning("Empty file: ", basename(f))
    return(NULL)
  }
  
  log2fc_col <- grep("DESeq2_log2FC", colnames(df), value = TRUE)
  
  if (length(log2fc_col) == 0) {
    warning("No DESeq2_log2FC column in ", basename(f))
    return(NULL)
  }
  
  df %>%
    select(gene_id, all_of(log2fc_col)) %>%
    rename(!!cond := all_of(log2fc_col))
})

# Remove NULL results
df_list <- df_list[!sapply(df_list, is.null)]

if (length(df_list) == 0) {
  stop("ERROR: No valid TSV files could be read")
}

# Merge and prepare matrix (LÓGICA ORIGINAL)
mat <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), df_list) %>%
  drop_na()

if (nrow(mat) == 0) {
  stop("ERROR: No common genes found across conditions for heatmap")
}

message("Heatmap will include ", nrow(mat), " genes across ", ncol(mat) - 1, " conditions")

mat_heat <- as.matrix(mat[, -1])
rownames(mat_heat) <- mat$gene_id
mat_heat_z <- t(scale(t(mat_heat)))

message("Generating integrated genes Heatmap...")

ph <- tryCatch({
  pheatmap(mat_heat_z,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           fontsize_row = 6)
}, error = function(e) {
  warning("Failed to generate heatmap: ", e$message)
  return(NULL)
})

if (!is.null(ph)) {
  outfile <- file.path(PLOTS_DIR, "Integrated_heatmap.png")
  
  tryCatch({
    ggsave(outfile, ph, width = 7, height = 6)
    message("Integrated genes Heatmap saved: ", outfile)
  }, error = function(e) {
    warning("Failed to save heatmap: ", e$message)
  })
}

# ==========================
# Integrated UpSet Plot
# ==========================
message("\n[4/4] === Integrated UpSet Plot ===")

# Reuse tsv_files already read (LÓGICA ORIGINAL)
gene_sets <- lapply(tsv_files, function(f) {
  cond <- sub("Integrated_(.*)_filtered\\.tsv", "\\1", basename(f))
  
  df <- tryCatch({
    read_tsv(f, show_col_types = FALSE)
  }, error = function(e) {
    warning("Failed to read ", basename(f), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  df$gene_id
})

names(gene_sets) <- sapply(tsv_files, function(f) {
  sub("Integrated_(.*)_filtered\\.tsv", "\\1", basename(f))
})

# Remove NULL results
gene_sets <- gene_sets[!sapply(gene_sets, is.null)]

if (length(gene_sets) == 0) {
  stop("ERROR: No valid gene sets for UpSet plot")
}

message("Generating UpSet plot with ", length(gene_sets), " conditions...")

outfile <- file.path(PLOTS_DIR, "Integrated_UpSet.png")

# VOLVER AL CÓDIGO ORIGINAL (sin tryCatch)
png(outfile, width = 1800, height = 1200, res = 150)

upset(
  fromList(gene_sets),
  nsets = length(gene_sets),
  order.by = "freq",
  sets.x.label = "Genes per condition",
  mainbar.y.label = "Intersection size"
)

dev.off()

message("UpSet plot saved: ", outfile)

# --------------------------
# Summary
# --------------------------
message("\n=== Visualization Summary ===")
message("PCA plots: 2 (IP, IN)")
message("Heatmap genes: ", nrow(mat))
message("UpSet conditions: ", length(gene_sets))
message("Output directory: ", PLOTS_DIR)

message("\n=== Visualization completed successfully ===")

sink('session_info.txt')
sessionInfo()
sink()
