# ==========================
# RIP-seq counts preparation
# Author: Adil Hannaoui Anaaoui
# ==========================

library(dplyr)
library(readr)

# --------------------------
# Load configuration
# --------------------------
source("R/config.R")

message("=== [Counts Preparation] Starting ===")

# --------------------------
# Validate required variables
# --------------------------
if (!exists("FEATURECOUNTS_DIR") || is.null(FEATURECOUNTS_DIR)) {
  stop("ERROR: FEATURECOUNTS_DIR not defined in config.R")
}

if (!dir.exists(FEATURECOUNTS_DIR)) {
  stop("ERROR: featureCounts directory not found at: ", FEATURECOUNTS_DIR)
}

if (!exists("CONDITIONS") || is.null(CONDITIONS) || length(CONDITIONS) == 0) {
  stop("ERROR: CONDITIONS not defined or empty in config.R")
}

if (!exists("DESEQ2_FC_DIR") || is.null(DESEQ2_FC_DIR)) {
  stop("ERROR: DESEQ2_FC_DIR not defined in config.R")
}

message("Using featureCounts directory: ", FEATURECOUNTS_DIR)
message("Expected conditions: ", paste(CONDITIONS, collapse = ", "))

# --------------------------
# 1. List all featureCounts files
# --------------------------
message("\n[1/6] Listing featureCounts files...")

files <- list.files(
  path = FEATURECOUNTS_DIR,
  pattern = "Counts_.*\\.txt$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("ERROR: No featureCounts files (Counts_*.txt) found in ", FEATURECOUNTS_DIR)
}

files <- sort(files)
message("Found ", length(files), " count files")

# --------------------------
# 2. Read all files into a list
# --------------------------
message("\n[2/6] Reading count files...")

counts_list <- list()
failed_reads <- character()

for (f in files) {
  tryCatch({
    # Read file
    df <- read_table(f, show_col_types = FALSE, col_types = "ci")
    
    # Validate columns
    if (ncol(df) < 2) {
      warning("File ", basename(f), " has fewer than 2 columns, skipping")
      failed_reads <- c(failed_reads, basename(f))
      next
    }
    
    # Ensure column names
    colnames(df)[1:2] <- c("Geneid", "counts")
    
    # Select and validate
    df <- df %>%
      dplyr::select(Geneid, counts) %>%
      filter(!is.na(Geneid), !is.na(counts))
    
    if (nrow(df) == 0) {
      warning("File ", basename(f), " is empty after filtering, skipping")
      failed_reads <- c(failed_reads, basename(f))
      next
    }
    
    # Extract sample name
    sample_name <- tools::file_path_sans_ext(basename(f))
    sample_name <- sub("^Counts_", "", sample_name)
    
    counts_list[[sample_name]] <- df
    message("  Read ", sample_name, ": ", nrow(df), " genes")
    
  }, error = function(e) {
    warning("Failed to read ", basename(f), ": ", e$message)
    failed_reads <<- c(failed_reads, basename(f))
  })
}

if (length(counts_list) == 0) {
  stop("ERROR: No count files could be read successfully")
}

if (length(failed_reads) > 0) {
  warning("Failed to read ", length(failed_reads), " files: ", 
          paste(failed_reads, collapse = ", "))
}

message("Successfully read ", length(counts_list), " count files")

# --------------------------
# 3. Identify common genes
# --------------------------
message("\n[3/6] Identifying common genes across samples...")

all_genes <- lapply(counts_list, `[[`, "Geneid")
common_genes <- Reduce(intersect, all_genes)

if (length(common_genes) == 0) {
  stop("ERROR: No common genes found across all samples")
}

# Report gene overlap statistics
total_genes <- unique(unlist(all_genes))
message("Total unique genes across all samples: ", length(total_genes))
message("Common genes across all samples: ", length(common_genes))
message("Percentage of genes in common: ", 
        round(100 * length(common_genes) / length(total_genes), 2), "%")

# Report per-sample gene counts
message("\nGenes per sample:")
for (sample in names(counts_list)) {
  n_genes <- nrow(counts_list[[sample]])
  n_common <- sum(counts_list[[sample]]$Geneid %in% common_genes)
  message("  ", sample, ": ", n_genes, " genes (", n_common, " common)")
}

# --------------------------
# 4. Build counts matrix
# --------------------------
message("\n[4/6] Building counts matrix...")

# Sort common genes for consistency
common_genes <- sort(common_genes)

counts_matrix <- sapply(counts_list, function(df) {
  df %>%
    filter(Geneid %in% common_genes) %>%
    arrange(Geneid) %>%  # Ensure same order
    pull(counts)
})

rownames(counts_matrix) <- common_genes

# Validate matrix
if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
  stop("ERROR: Empty counts matrix generated")
}

if (any(is.na(counts_matrix))) {
  warning("Counts matrix contains ", sum(is.na(counts_matrix)), " NA values")
  # Replace NAs with 0
  counts_matrix[is.na(counts_matrix)] <- 0
}

message("Counts matrix dimensions: ", nrow(counts_matrix), " genes x ", 
        ncol(counts_matrix), " samples")

# Basic statistics
message("\nCounts matrix statistics:")
message("  Min count: ", min(counts_matrix))
message("  Max count: ", max(counts_matrix))
message("  Median count: ", median(counts_matrix))
message("  Total counts: ", sum(counts_matrix))

# --------------------------
# 5. Define sample metadata
# --------------------------
message("\n[5/6] Creating sample metadata...")

sample_names <- colnames(counts_matrix)

# Auto-detect conditions from sample names if not matching
if (length(sample_names) != length(CONDITIONS)) {
  warning("Number of samples (", length(sample_names), 
          ") does not match number of CONDITIONS (", length(CONDITIONS), ")")
  warning("Attempting to auto-detect conditions from sample names...")
  
  # Extract condition from sample name pattern
  detected_conditions <- character(length(sample_names))
  
  for (i in seq_along(sample_names)) {
    sample <- sample_names[i]
    
    # Try to match each known condition
    matched <- FALSE
    for (cond in CONDITIONS) {
      if (grepl(cond, sample, ignore.case = TRUE)) {
        detected_conditions[i] <- cond
        matched <- TRUE
        break
      }
    }
    
    if (!matched) {
      # Extract prefix before _IP or _IN
      detected_conditions[i] <- sub("_[Ii][PpNn].*", "", sample)
    }
  }
  
  message("Auto-detected conditions:")
  print(table(detected_conditions))
  
  colData <- data.frame(
    sample = sample_names,
    condition = factor(detected_conditions)
  )
  
} else {
  # Original behavior: use CONDITIONS directly
  colData <- data.frame(
    sample = sample_names,
    condition = factor(CONDITIONS)
  )
}

rownames(colData) <- sample_names

# Validate colData
if (nrow(colData) == 0) {
  stop("ERROR: Empty sample metadata generated")
}

if (any(is.na(colData$condition))) {
  stop("ERROR: Some samples have NA condition assignment")
}

message("\nSample metadata summary:")
print(table(colData$condition))

message("\nSample-condition mapping:")
print(colData)

# --------------------------
# 6. Save outputs
# --------------------------
message("\n[6/6] Saving processed data...")

dir.create(DESEQ2_FC_DIR, showWarnings = FALSE, recursive = TRUE)

# Save counts matrix
counts_file <- file.path(DESEQ2_FC_DIR, "counts_matrix.rds")
saveRDS(counts_matrix, file = counts_file)
message("  Saved counts matrix: ", counts_file)

# Save metadata
coldata_file <- file.path(DESEQ2_FC_DIR, "colData.rds")
saveRDS(colData, file = coldata_file)
message("  Saved sample metadata: ", coldata_file)

# Also save as TSV for easy inspection
counts_tsv <- file.path(DESEQ2_FC_DIR, "counts_matrix.tsv")
write_tsv(
  as.data.frame(counts_matrix) %>% 
    tibble::rownames_to_column("Geneid"),
  counts_tsv
)
message("  Saved counts matrix (TSV): ", counts_tsv)

coldata_tsv <- file.path(DESEQ2_FC_DIR, "colData.tsv")
# FIX: colData ya tiene la columna "sample", así que no usar rownames_to_column
# En su lugar, simplemente guardar colData directamente
write_tsv(colData, coldata_tsv)
message("  Saved sample metadata (TSV): ", coldata_tsv)

# --------------------------
# Summary report
# --------------------------
message("\n=== Counts Preparation Summary ===")
message("Samples processed: ", ncol(counts_matrix))
message("Genes retained: ", nrow(counts_matrix))
message("Conditions: ", paste(unique(colData$condition), collapse = ", "))
message("Output directory: ", DESEQ2_FC_DIR)
message("\n=== Counts preparation completed successfully ===")
