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

# --------------------------
# 1. List all featureCounts files
# --------------------------
message("Listing featureCounts files...")

files <- list.files(
  path = FEATURECOUNTS_DIR,
  pattern = "Counts_.*\\.txt$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No featureCounts files found in ", FEATURECOUNTS_DIR)
}

files <- sort(files)

# --------------------------
# 2. Read all files into a list
# --------------------------
message("Reading count files...")

counts_list <- lapply(files, function(f) {
  read_table(f, show_col_types = FALSE) %>%
    dplyr::select(Geneid, counts = 2)
})

names(counts_list) <- basename(files) %>%
  tools::file_path_sans_ext()

# --------------------------
# 3. Identify common genes
# --------------------------
common_genes <- Reduce(intersect, lapply(counts_list, `[[`, "Geneid"))

if (length(common_genes) == 0) {
  stop("No common genes found across samples")
}

# --------------------------
# 4. Build counts matrix
# --------------------------
counts_matrix <- sapply(counts_list, function(df) {
  df %>%
    filter(Geneid %in% common_genes) %>%
    pull(counts)
})

rownames(counts_matrix) <- common_genes

# --------------------------
# 5. Define sample metadata
# --------------------------
sample_names <- colnames(counts_matrix)

if (length(sample_names) != length(CONDITIONS)) {
  stop("Number of samples does not match number of conditions")
}

colData <- data.frame(
  sample = sample_names,
  condition = factor(CONDITIONS)
)
table(colData$sample, colData$condition)
rownames(colData) <- sample_names

# --------------------------
# 6. Save outputs
# --------------------------
message("Saving processed data...")

saveRDS(counts_matrix, file = file.path(OUTPUT_DIR, "counts_matrix.rds"))
saveRDS(colData, file = file.path(OUTPUT_DIR, "colData.rds"))

message("Counts preparation completed successfully")
