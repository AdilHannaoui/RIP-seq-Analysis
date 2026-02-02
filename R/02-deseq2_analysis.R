# ==========================
# RIP-seq DESeq2 analysis
# Author: Adil Hannaoui Anaaoui
# ==========================

# --------------------------
# Load configuration
# --------------------------
source("R/config.R") 

library(DESeq2)
library(dplyr)

# --------------------------
# Load data
# --------------------------
counts_matrix <- readRDS(COUNTS_MATRIX_PATH)    # counts matrix: genes x samples
colData <- readRDS(SAMPLE_METADATA_PATH)       # metadata: samples x conditions

# --------------------------
# Create DESeqDataSet
# --------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~ condition
)

# --------------------------
# Filter low-count genes
# --------------------------
dds <- dds[rowSums(counts(dds)) > MIN_COUNTS_FILTER, ]

# --------------------------
# Set reference level
# --------------------------
dds$condition <- relevel(dds$condition, ref = REFERENCE_CONDITION)

# --------------------------
# Run DESeq2
# --------------------------
dds <- DESeq(dds)

# --------------------------
# Extract results for contrasts
# --------------------------
res_WT_vs_M1 <- results(dds, contrast = c("condition", "M1", REFERENCE_CONDITION))
res_WT_vs_M12 <- results(dds, contrast = c("condition", "M12", REFERENCE_CONDITION))

# --------------------------
# Filter significant genes
# --------------------------
res_WT_vs_M1_sig <- res_WT_vs_M1[which(res_WT_vs_M1$padj < PADJ_THRESHOLD), ]
res_WT_vs_M12_sig <- res_WT_vs_M12[which(res_WT_vs_M12$padj < PADJ_THRESHOLD), ]

# --------------------------
# Save results
# --------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(dds, file = file.path(OUTPUT_DIR, "dds.rds"))
saveRDS(res_WT_vs_M1_sig, file = file.path(OUTPUT_DIR, "DESeq2_WT_vs_M1_sig.rds"))
saveRDS(res_WT_vs_M12_sig, file = file.path(OUTPUT_DIR, "DESeq2_WT_vs_M12_sig.rds"))

write.csv(as.data.frame(res_WT_vs_M1_sig), file = file.path(OUTPUT_DIR, "DESeq2_WT_vs_M1_sig.csv"))
write.csv(as.data.frame(res_WT_vs_M12_sig), file = file.path(OUTPUT_DIR, "DESeq2_WT_vs_M12_sig.csv"))

cat("DESeq2 analysis completed. Significant results saved in:", OUTPUT_DIR, "\n")
