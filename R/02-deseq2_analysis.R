# ==========================
# RIP-seq DESeq2 analysis (featureCounts)
# Author: Adil Hannaoui Anaaoui
# ==========================

library(DESeq2)
library(dplyr)

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

message("=== [DESeq2 FeatureCounts] Starting ===")

# --------------------------
# Validate required variables
# --------------------------
if (!exists("COUNTS_MATRIX_PATH_FC") || is.null(COUNTS_MATRIX_PATH_FC)) {
  stop("ERROR: COUNTS_MATRIX_PATH_FC not defined in config.R")
}

if (!exists("SAMPLE_METADATA_PATH_FC") || is.null(SAMPLE_METADATA_PATH_FC)) {
  stop("ERROR: SAMPLE_METADATA_PATH_FC not defined in config.R")
}

if (!exists("MIN_COUNTS_FILTER") || is.null(MIN_COUNTS_FILTER)) {
  MIN_COUNTS_FILTER <- 10
  warning("MIN_COUNTS_FILTER not defined, using default: ", MIN_COUNTS_FILTER)
}

if (!exists("PADJ_THRESHOLD") || is.null(PADJ_THRESHOLD)) {
  PADJ_THRESHOLD <- 0.05
  warning("PADJ_THRESHOLD not defined, using default: ", PADJ_THRESHOLD)
}

if (!exists("DESEQ2_FC_DIR") || is.null(DESEQ2_FC_DIR)) {
  stop("ERROR: DESEQ2_FC_DIR not defined in config.R")
}

message("Using parameters:")
message("  MIN_COUNTS_FILTER: ", MIN_COUNTS_FILTER)
message("  PADJ_THRESHOLD: ", PADJ_THRESHOLD)

# --------------------------
# 1. Load data
# --------------------------
message("\n[1/6] Loading counts matrix and metadata...")

if (!file.exists(COUNTS_MATRIX_PATH_FC)) {
  stop("ERROR: Counts matrix not found at: ", COUNTS_MATRIX_PATH_FC)
}

if (!file.exists(SAMPLE_METADATA_PATH_FC)) {
  stop("ERROR: Sample metadata not found at: ", SAMPLE_METADATA_PATH_FC)
}

counts_matrix <- readRDS(COUNTS_MATRIX_PATH_FC)
colData <- readRDS(SAMPLE_METADATA_PATH_FC)

# Validate data
if (!is.matrix(counts_matrix) && !is.data.frame(counts_matrix)) {
  stop("ERROR: counts_matrix is not a matrix or data frame")
}

if (!is.data.frame(colData)) {
  stop("ERROR: colData is not a data frame")
}

if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
  stop("ERROR: Empty counts matrix")
}

if (nrow(colData) == 0) {
  stop("ERROR: Empty sample metadata")
}

if (ncol(counts_matrix) != nrow(colData)) {
  stop("ERROR: Number of samples in counts matrix (", ncol(counts_matrix), 
       ") does not match metadata (", nrow(colData), ")")
}

message("Loaded counts matrix: ", nrow(counts_matrix), " genes x ", 
        ncol(counts_matrix), " samples")
message("Loaded metadata: ", nrow(colData), " samples")

# Check for 'condition' column
if (!"condition" %in% colnames(colData)) {
  stop("ERROR: 'condition' column not found in metadata")
}

message("\nConditions in metadata:")
print(table(colData$condition))

# --------------------------
# 2. Create DESeqDataSet
# --------------------------
message("\n[2/6] Creating DESeqDataSet...")

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~ condition
)

message("Initial dataset: ", nrow(dds), " genes, ", ncol(dds), " samples")

# --------------------------
# 3. Filter low-count genes
# --------------------------
message("\n[3/6] Filtering low-count genes...")

genes_before <- nrow(dds)
dds <- dds[rowSums(counts(dds)) > MIN_COUNTS_FILTER, ]
genes_after <- nrow(dds)

message("Removed ", genes_before - genes_after, " genes with rowSums <= ", 
        MIN_COUNTS_FILTER)
message("Retained ", genes_after, " genes for analysis")

if (nrow(dds) == 0) {
  stop("ERROR: No genes remain after filtering")
}

# --------------------------
# 4. Run DESeq2
# --------------------------
message("\n[4/6] Running DESeq2 differential expression analysis...")

dds <- DESeq(dds)

message("DESeq2 analysis completed")

# --------------------------
# 5. Extract RIP-seq contrasts
# --------------------------
message("\n[5/6] Extracting differential expression results...")

# Get available conditions
available_conditions <- levels(colData$condition)
message("Available conditions: ", paste(available_conditions, collapse = ", "))

# Auto-detect IP/IN pairs
condition_pairs <- list()

# Extract unique sample groups (M1, M12, WT, etc.)
sample_groups <- unique(gsub("_[Ii][PpNn].*", "", available_conditions))

message("Detected sample groups: ", paste(sample_groups, collapse = ", "))

# Store results
results_list <- list()
sig_results_list <- list()

for (group in sample_groups) {
  # Find IP and IN conditions for this group (case-insensitive)
  ip_cond <- grep(paste0(group, "_IP"), available_conditions, 
                  ignore.case = TRUE, value = TRUE)
  in_cond <- grep(paste0(group, "_IN"), available_conditions, 
                  ignore.case = TRUE, value = TRUE)
  
  if (length(ip_cond) == 0 || length(in_cond) == 0) {
    warning("Skipping ", group, ": missing IP or IN condition")
    next
  }
  
  # Use first match if multiple
  if (length(ip_cond) > 1) {
    warning("Multiple IP conditions found for ", group, ", using first: ", ip_cond[1])
    ip_cond <- ip_cond[1]
  }
  
  if (length(in_cond) > 1) {
    warning("Multiple IN conditions found for ", group, ", using first: ", in_cond[1])
    in_cond <- in_cond[1]
  }
  
  message("  Extracting results for ", group, " (", ip_cond, " vs ", in_cond, ")...")
  
  # Extract results
  tryCatch({
    res <- results(dds, contrast = c("condition", ip_cond, in_cond))
    
    # Convert to data frame and filter significant
    res_df <- as.data.frame(res)
    res_sig <- res_df %>% filter(!is.na(padj) & padj < PADJ_THRESHOLD)
    
    # Store results
    results_list[[group]] <- res_df
    sig_results_list[[group]] <- res_sig
    
    message("    Total genes tested: ", sum(!is.na(res_df$padj)))
    message("    Significant genes (padj < ", PADJ_THRESHOLD, "): ", nrow(res_sig))
    
    if (nrow(res_sig) > 0) {
      message("    Upregulated (log2FC > 0): ", sum(res_sig$log2FoldChange > 0))
      message("    Downregulated (log2FC < 0): ", sum(res_sig$log2FoldChange < 0))
    }
    
  }, error = function(e) {
    warning("Failed to extract results for ", group, ": ", e$message)
  })
}

if (length(results_list) == 0) {
  stop("ERROR: No results could be extracted for any sample group")
}

if (length(sig_results_list) == 0) {
  warning("No significant results found for any sample group")
}

# --------------------------
# 6. Save results
# --------------------------
message("\n[6/6] Saving results...")

dir.create(DESEQ2_FC_DIR, showWarnings = FALSE, recursive = TRUE)

# Save DESeq2 object
dds_file <- file.path(DESEQ2_FC_DIR, "dds.rds")
saveRDS(dds, file = dds_file)
message("  Saved DESeq2 object: ", basename(dds_file))

# Save results for each group
for (group in names(results_list)) {
  # Full results
  full_rds <- file.path(DESEQ2_FC_DIR, 
                        paste0("DESeq2_", group, "_IP_vs_IN_full.rds"))
  saveRDS(results_list[[group]], file = full_rds)
  
  full_csv <- file.path(DESEQ2_FC_DIR, 
                        paste0("DESeq2_", group, "_IP_vs_IN_full.csv"))
  write.csv(results_list[[group]], file = full_csv)
  
  message("  Saved ", group, " full results: ", basename(full_rds), ", ", 
          basename(full_csv))
}

# Save significant results (if any)
for (group in names(sig_results_list)) {
  if (nrow(sig_results_list[[group]]) > 0) {
    sig_rds <- file.path(DESEQ2_FC_DIR, 
                         paste0("DESeq2_", group, "_IP_vs_IN_sig.rds"))
    saveRDS(sig_results_list[[group]], file = sig_rds)
    
    sig_csv <- file.path(DESEQ2_FC_DIR, 
                         paste0("DESeq2_", group, "_IP_vs_IN_sig.csv"))
    write.csv(sig_results_list[[group]], file = sig_csv)
    
    message("  Saved ", group, " significant results: ", basename(sig_rds), ", ", 
            basename(sig_csv))
  }
}

# --------------------------
# Summary report
# --------------------------
message("\n=== DESeq2 FeatureCounts Analysis Summary ===")
message("Sample groups analyzed: ", paste(names(results_list), collapse = ", "))
message("Genes in final dataset: ", nrow(dds))
message("Output directory: ", DESEQ2_FC_DIR)

if (length(sig_results_list) > 0) {
  message("\nSignificant genes per group (padj < ", PADJ_THRESHOLD, "):")
  for (group in names(sig_results_list)) {
    n_sig <- nrow(sig_results_list[[group]])
    message("  ", group, ": ", n_sig, " genes")
  }
} else {
  message("\nNo significant genes found in any group (padj < ", PADJ_THRESHOLD, ")")
}

message("\n=== DESeq2 analysis completed successfully ===")