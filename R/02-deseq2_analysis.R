# ==========================
# RIP-seq DESeq2 analysis
# Author: Adil Hannaoui Anaaoui
# ==========================

source("R/config.R")

library(DESeq2)
library(dplyr)

# --------------------------
# Load data
# --------------------------
counts_matrix <- readRDS(COUNTS_MATRIX_PATH)
colData <- readRDS(SAMPLE_METADATA_PATH)

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
# Run DESeq2
# --------------------------
dds <- DESeq(dds)

# --------------------------
# Extract RIP-seq contrasts
# --------------------------
res_M1_IP_vs_IN   <- results(dds, contrast = c("condition", "M1_IP", "M1_IN"))
res_M12_IP_vs_IN  <- results(dds, contrast = c("condition", "M12_IP", "M12_IN"))
res_WT_IP_vs_IN   <- results(dds, contrast = c("condition", "WT_IP", "WT_IN"))

# --------------------------
# Filter significant genes
# --------------------------
sig_M1  <- res_M1_IP_vs_IN  %>% as.data.frame() %>% filter(padj < PADJ_THRESHOLD)
sig_M12 <- res_M12_IP_vs_IN %>% as.data.frame() %>% filter(padj < PADJ_THRESHOLD)
sig_WT  <- res_WT_IP_vs_IN  %>% as.data.frame() %>% filter(padj < PADJ_THRESHOLD)

# --------------------------
# Save results
# --------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(dds, file = file.path(OUTPUT_DIR, "dds.rds"))

saveRDS(sig_M1,  file = file.path(OUTPUT_DIR, "DESeq2_M1_IP_vs_IN_sig.rds"))
saveRDS(sig_M12, file = file.path(OUTPUT_DIR, "DESeq2_M12_IP_vs_IN_sig.rds"))
saveRDS(sig_WT,  file = file.path(OUTPUT_DIR, "DESeq2_WT_IP_vs_IN_sig.rds"))

write.csv(sig_M1,  file = file.path(OUTPUT_DIR, "DESeq2_M1_IP_vs_IN_sig.csv"))
write.csv(sig_M12, file = file.path(OUTPUT_DIR, "DESeq2_M12_IP_vs_IN_sig.csv"))
write.csv(sig_WT,  file = file.path(OUTPUT_DIR, "DESeq2_WT_IP_vs_IN_sig.csv"))

cat("DESeq2 analysis completed. Results saved in:", OUTPUT_DIR, "\n")

write.csv(as.data.frame(res_WT_vs_M1_sig), file = file.path(OUTPUT_DIR, "DESeq2_WT_vs_M1_sig.csv"))
write.csv(as.data.frame(res_WT_vs_M12_sig), file = file.path(OUTPUT_DIR, "DESeq2_WT_vs_M12_sig.csv"))

cat("DESeq2 analysis completed. Significant results saved in:", OUTPUT_DIR, "\n")
