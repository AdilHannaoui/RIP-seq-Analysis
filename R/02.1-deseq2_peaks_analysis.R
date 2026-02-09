# ==========================
# RIP-seq DESeq2 analysis (MACS3 peaks)
# Author: Adil Hannaoui Anaaoui
# ==========================

library(DESeq2)
library(dplyr)
library(readr)

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

message("=== [DESeq2 MACS3] Starting ===")

# --------------------------
# Validate required variables
# --------------------------
if (!exists("MIN_COUNTS_FILTER") || is.null(MIN_COUNTS_FILTER)) {
  MIN_COUNTS_FILTER <- 10
  warning("MIN_COUNTS_FILTER not defined, using default: ", MIN_COUNTS_FILTER)
}

if (!exists("PADJ_THRESHOLD") || is.null(PADJ_THRESHOLD)) {
  PADJ_THRESHOLD <- 0.05
  warning("PADJ_THRESHOLD not defined, using default: ", PADJ_THRESHOLD)
}

if (!exists("DESEQ2_MACS3_DIR") || is.null(DESEQ2_MACS3_DIR)) {
  stop("ERROR: DESEQ2_MACS3_DIR not defined in config.R")
}

message("Using parameters:")
message("  MIN_COUNTS_FILTER: ", MIN_COUNTS_FILTER)
message("  PADJ_THRESHOLD: ", PADJ_THRESHOLD)

# --------------------------
# 1. Load MACS3 common counts files
# --------------------------
message("\n[1/7] Loading MACS3 common counts files...")

MACS3_DIR <- "output/macs3"

if (!dir.exists(MACS3_DIR)) {
  stop("ERROR: MACS3 directory not found at: ", MACS3_DIR)
}

COUNTS_FILES <- list.files(
  path = MACS3_DIR,
  pattern = "_common_counts\\.txt$",
  full.names = TRUE
)

if (length(COUNTS_FILES) == 0) {
  stop("ERROR: No common counts files (*_common_counts.txt) found in ", MACS3_DIR)
}

message("Found ", length(COUNTS_FILES), " count files:")
for (f in COUNTS_FILES) {
  message("  - ", basename(f))
}

# --------------------------
# 2. Read and process count files
# --------------------------
message("\n[2/7] Reading and processing count files...")

count_list <- list()
failed_reads <- character()

for (f in COUNTS_FILES) {
  tryCatch({
    df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Validate expected columns
    expected_cols <- c("chrom", "start", "end", "name", "score", "strand", "Gene_id", "Biotype")
    if (!all(expected_cols %in% colnames(df))) {
      warning("File ", basename(f), " missing expected columns, skipping")
      failed_reads <- c(failed_reads, basename(f))
      next
    }
    
    # Extract sample name
    sample_name <- sub("_common_counts\\.txt$", "", basename(f))
    
    # Identify count columns (typically IP1, IN1, IP2, IN2, ...)
    count_cols <- setdiff(colnames(df), expected_cols)
    
    if (length(count_cols) == 0) {
      warning("File ", basename(f), " has no count columns, skipping")
      failed_reads <- c(failed_reads, basename(f))
      next
    }
    
    # Rename count columns to include sample name
    new_count_names <- paste0(sample_name, "_", count_cols)
    colnames(df)[colnames(df) %in% count_cols] <- new_count_names
    
    count_list[[sample_name]] <- df
    message("  Read ", sample_name, ": ", nrow(df), " peaks, ", length(count_cols), " count columns")
    
  }, error = function(e) {
    warning("Failed to read ", basename(f), ": ", e$message)
    failed_reads <<- c(failed_reads, basename(f))
  })
}

if (length(count_list) == 0) {
  stop("ERROR: No count files could be read successfully")
}

if (length(failed_reads) > 0) {
  warning("Failed to read ", length(failed_reads), " files: ", 
          paste(failed_reads, collapse = ", "))
}

message("Successfully processed ", length(count_list), " count files")

# --------------------------
# 3. Merge all count tables
# --------------------------
message("\n[3/7] Merging count tables...")

merge_cols <- c("chrom", "start", "end", "name", "score", "strand", "Gene_id", "Biotype")

counts_merged <- Reduce(function(x, y) {
  merge(x, y, by = merge_cols, all = TRUE)
}, count_list)

if (is.null(counts_merged) || nrow(counts_merged) == 0) {
  stop("ERROR: Merged counts table is empty")
}

message("Merged table: ", nrow(counts_merged), " peaks")

# Replace NAs with 0 (peaks not present in all samples)
count_cols <- setdiff(colnames(counts_merged), merge_cols)
counts_merged[, count_cols][is.na(counts_merged[, count_cols])] <- 0

na_count <- sum(is.na(counts_merged[, count_cols]))
if (na_count > 0) {
  message("Replaced ", na_count, " NA values with 0")
}

# --------------------------
# 4. Collapse peaks by gene (SUM)
# --------------------------
message("\n[4/7] Collapsing peaks by gene...")

count_data <- counts_merged %>%
  select(Gene_id, all_of(count_cols)) %>%
  group_by(Gene_id) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

# Convert to matrix
counts_matrix <- as.data.frame(count_data)
rownames(counts_matrix) <- counts_matrix$Gene_id
counts_matrix$Gene_id <- NULL

# Validate matrix
if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
  stop("ERROR: Empty counts matrix after collapsing by gene")
}

message("Counts matrix: ", nrow(counts_matrix), " genes x ", ncol(counts_matrix), " samples")

# Remove genes with all zeros
all_zero_genes <- rowSums(counts_matrix) == 0
if (sum(all_zero_genes) > 0) {
  message("Removing ", sum(all_zero_genes), " genes with all zero counts")
  counts_matrix <- counts_matrix[!all_zero_genes, ]
}

# --------------------------
# 5. Build colData automatically
# --------------------------
message("\n[5/7] Creating sample metadata...")

sample_names <- colnames(counts_matrix)

# Extract sample group (M1, M12, WT) and condition (IP/IN)
sample_group <- sub("_.*", "", sample_names)
condition <- ifelse(grepl("_IP", sample_names), "IP", "IN")

# Auto-detect available groups
available_groups <- unique(sample_group)
message("Detected sample groups: ", paste(available_groups, collapse = ", "))

colData <- data.frame(
  row.names = sample_names,
  sample_group = factor(sample_group, levels = sort(available_groups)),
  condition = factor(condition, levels = c("IN", "IP"))
)

# Create combined factor for contrast
colData$group_condition <- factor(paste0(colData$sample_group, "_", colData$condition))

message("\nSample metadata summary:")
print(table(colData$sample_group, colData$condition))

# --------------------------
# 6. Create DESeqDataSet and run DESeq2
# --------------------------
message("\n[6/7] Running DESeq2 analysis...")

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~ group_condition
)

message("Initial dataset: ", nrow(dds), " genes")

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > MIN_COUNTS_FILTER, ]
message("After filtering (rowSums > ", MIN_COUNTS_FILTER, "): ", nrow(dds), " genes")

if (nrow(dds) == 0) {
  stop("ERROR: No genes remain after filtering")
}

# Run DESeq2
message("Running DESeq2 differential expression analysis...")
dds <- DESeq(dds)

# --------------------------
# 7. Extract results: IP vs IN per sample group
# --------------------------
message("\n[7/7] Extracting differential expression results...")

# Create list to store results
results_list <- list()
sig_results_list <- list()

# Extract results for each available group
for (group in available_groups) {
  ip_name <- paste0(group, "_IP")
  in_name <- paste0(group, "_IN")
  
  # Check if both conditions exist
  if (!ip_name %in% levels(colData$group_condition) || 
      !in_name %in% levels(colData$group_condition)) {
    warning("Skipping ", group, ": missing IP or IN samples")
    next
  }
  
  message("  Extracting results for ", group, " (IP vs IN)...")
  
  res <- results(dds, contrast = c("group_condition", ip_name, in_name))
  
  # Filter significant results
  res_sig <- res[which(!is.na(res$padj) & res$padj < PADJ_THRESHOLD), ]
  
  # Store results
  results_list[[group]] <- res
  sig_results_list[[group]] <- res_sig
  
  message("    Total genes tested: ", sum(!is.na(res$padj)))
  message("    Significant genes (padj < ", PADJ_THRESHOLD, "): ", nrow(res_sig))
  
  if (nrow(res_sig) > 0) {
    message("    Upregulated (log2FC > 0): ", sum(res_sig$log2FoldChange > 0))
    message("    Downregulated (log2FC < 0): ", sum(res_sig$log2FoldChange < 0))
  }
}

if (length(sig_results_list) == 0) {
  warning("No significant results found for any sample group")
}

# --------------------------
# 8. Save results
# --------------------------
message("\nSaving results...")

dir.create(DESEQ2_MACS3_DIR, showWarnings = FALSE, recursive = TRUE)

# Core objects for downstream modules
saveRDS(dds, file = file.path(DESEQ2_MACS3_DIR, "dds_macs3.rds"))
message("  Saved DESeq2 object: dds_macs3.rds")

saveRDS(colData, file = file.path(DESEQ2_MACS3_DIR, "colData_macs3.rds"))
message("  Saved metadata: colData_macs3.rds")

saveRDS(counts_matrix, file = file.path(DESEQ2_MACS3_DIR, "counts_matrix_macs3.rds"))
message("  Saved counts matrix: counts_matrix_macs3.rds")

# Save results for each group
for (group in names(sig_results_list)) {
  # RDS files
  rds_file <- file.path(DESEQ2_MACS3_DIR, paste0("DESeq2_res_", group, "_sig.rds"))
  saveRDS(sig_results_list[[group]], file = rds_file)
  
  # CSV files
  csv_file <- file.path(DESEQ2_MACS3_DIR, paste0("DESeq2_res_", group, "_sig.csv"))
  write.csv(as.data.frame(sig_results_list[[group]]), file = csv_file)
  
  message("  Saved ", group, " results: ", basename(rds_file), ", ", basename(csv_file))
}

# Also save full (non-filtered) results
for (group in names(results_list)) {
  csv_file <- file.path(DESEQ2_MACS3_DIR, paste0("DESeq2_res_", group, "_full.csv"))
  write.csv(as.data.frame(results_list[[group]]), file = csv_file)
  message("  Saved ", group, " full results: ", basename(csv_file))
}

# --------------------------
# Summary report
# --------------------------
message("\n=== DESeq2 MACS3 Analysis Summary ===")
message("Sample groups analyzed: ", paste(names(results_list), collapse = ", "))
message("Genes in final dataset: ", nrow(dds))
message("Output directory: ", DESEQ2_MACS3_DIR)

message("\nSignificant genes per group (padj < ", PADJ_THRESHOLD, "):")
for (group in names(sig_results_list)) {
  n_sig <- nrow(sig_results_list[[group]])
  message("  ", group, ": ", n_sig, " genes")
}

message("\n=== DESeq2 analysis completed successfully ===")