############################################################
# RNA-seq Differential Expression: GSE81089 Tumor vs Normal
# Reconstructed from the user's R history.
#
# IMPORTANT:
# - Uses the raw featureCounts file, NOT the FPKM file.
# - Checks sample names and paired-sample metadata before DESeq2.
# - Review printed checks before interpreting results.
# - Outputs are written to an "results_rebuilt" subfolder.
############################################################

# 0. Packages ------------------------------------------------
required_cran <- c("ggplot2", "pheatmap")
required_bioc <- c(
  "DESeq2", "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "enrichplot"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(pheatmap)
})

# 1. Input and output paths ----------------------------------
# Change this path if the dataset is stored elsewhere.
input_file <- "C:/Users/Lenovo/Downloads/GSE81089_readcounts_featurecounts.tsv.gz"

if (!file.exists(input_file)) {
  stop(
    "Raw count file not found. Update 'input_file' to the correct path: ",
    input_file
  )
}

output_dir <- file.path(
  "C:/Users/Lenovo/OneDrive/Documents",
  "RNAseq_Tumor_vs_Normal",
  "results_rebuilt"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Read raw featureCounts data -----------------------------
counts_table <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Input dimensions:", nrow(counts_table), "rows x",
    ncol(counts_table), "columns\n")
print(head(counts_table[, seq_len(min(6, ncol(counts_table)))]))

# featureCounts commonly includes annotation columns before sample counts.
# This history indicates an Ensembl ID column named Ensembl_gene_id.
if (!("Ensembl_gene_id" %in% colnames(counts_table))) {
  stop(
    "Expected column 'Ensembl_gene_id' was not found. Inspect colnames(counts_table) ",
    "and adjust the gene-ID column/sample-column selection."
  )
}

gene_ids <- as.character(counts_table$Ensembl_gene_id)

# Remove featureCounts summary rows and empty IDs.
summary_rows <- gene_ids %in% c(
  "__no_feature", "__ambiguous", "__too_low_aQual",
  "__not_aligned", "__alignment_not_unique", ""
) | is.na(gene_ids)

# History suggests count columns begin after the gene ID and contain sample names.
# To avoid accidentally including feature annotation columns, select columns
# whose names look like the project's L...T / L...N sample IDs.
sample_cols <- grep("^L.*[NT]$", colnames(counts_table), value = TRUE)

if (length(sample_cols) == 0) {
  stop(
    "No sample columns matching '^L.*[NT]$' were found. Inspect colnames(counts_table) ",
    "and update sample_cols before proceeding."
  )
}

count_df <- counts_table[, sample_cols, drop = FALSE]
rownames(count_df) <- gene_ids

# Remove summary rows and invalid gene IDs.
valid_rows <- !summary_rows & !duplicated(gene_ids)
count_df <- count_df[valid_rows, , drop = FALSE]

# Convert counts to numeric and validate.
count_matrix <- as.matrix(
  data.frame(lapply(as.data.frame(count_df), function(x) as.numeric(as.character(x))),
             check.names = FALSE)
)
rownames(count_matrix) <- rownames(count_df)

if (anyNA(count_matrix)) {
  stop("NA values found in the count matrix. Inspect the source file before continuing.")
}
if (any(count_matrix < 0)) {
  stop("Negative values found. Confirm that the input is raw integer counts, not FPKM.")
}
if (any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
  stop("Non-integer values found. DESeq2 requires raw counts, not FPKM/normalized values.")
}
storage.mode(count_matrix) <- "integer"

cat("Count matrix dimensions after cleaning:", nrow(count_matrix), "genes x",
    ncol(count_matrix), "samples\n")
cat("Missing count values:", sum(is.na(count_matrix)), "\n")

# 3. Build sample metadata -----------------------------------
sample_names <- colnames(count_matrix)
condition <- ifelse(grepl("N$", sample_names), "Normal", "Tumor")

# Matches the patient-ID logic found in the saved history.
patient_id <- sub("[NT]$", "", sample_names)
patient_id[patient_id == "L608T_2122"] <- "L608"
patient_id[patient_id == "L771T_1"] <- "L771"

metadata <- data.frame(
  sample = sample_names,
  condition = factor(condition, levels = c("Normal", "Tumor")),
  patient_id = patient_id,
  row.names = sample_names,
  stringsAsFactors = FALSE
)

cat("\nSamples by condition:\n")
print(table(metadata$condition))
cat("\nPatient x condition table:\n")
print(table(metadata$patient_id, metadata$condition))

# Keep only patients represented in both conditions.
patient_table <- table(metadata$patient_id, metadata$condition)
paired_patients <- rownames(
  patient_table[
    patient_table[, "Normal"] > 0 & patient_table[, "Tumor"] > 0,
    ,
    drop = FALSE
  ]
)

if (length(paired_patients) == 0) {
  stop("No paired patients detected. Review sample names and metadata construction.")
}

paired_metadata <- metadata[metadata$patient_id %in% paired_patients, , drop = FALSE]
paired_metadata <- paired_metadata[
  order(paired_metadata$patient_id, paired_metadata$condition),
  ,
  drop = FALSE
]
count_matrix_paired <- count_matrix[, paired_metadata$sample, drop = FALSE]

stopifnot(identical(colnames(count_matrix_paired), paired_metadata$sample))

cat("\nPaired patients:", length(paired_patients), "\n")
cat("Paired samples by condition:\n")
print(table(paired_metadata$condition))

write.csv(metadata, file.path(output_dir, "sample_metadata_all.csv"), row.names = FALSE)
write.csv(paired_metadata, file.path(output_dir, "sample_metadata_paired.csv"), row.names = FALSE)

# 4. Filter low-count genes and run DESeq2 --------------------
keep <- rowSums(count_matrix_paired) >= 10
count_matrix_filtered <- count_matrix_paired[keep, , drop = FALSE]
cat("Genes before filtering:", nrow(count_matrix_paired), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = paired_metadata,
  design = ~ patient_id + condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)

write.csv(res_df, file.path(output_dir, "DESeq2_all_results.csv"), row.names = FALSE)

# 5. Significant DEGs ----------------------------------------
# Main threshold used in the history: adjusted p-value < 0.05 and |LFC| > 1.
sig <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) > 1,
  ,
  drop = FALSE
]

up_genes <- sig[sig$log2FoldChange > 0, , drop = FALSE]
down_genes <- sig[sig$log2FoldChange < 0, , drop = FALSE]

cat("Significant DEGs:", nrow(sig), "\n")
cat("Upregulated:", nrow(up_genes), "\n")
cat("Downregulated:", nrow(down_genes), "\n")

write.csv(sig, file.path(output_dir, "DESeq2_significant_genes_rebuilt.csv"), row.names = FALSE)
write.csv(up_genes, file.path(output_dir, "Upregulated_genes_rebuilt.csv"), row.names = FALSE)
write.csv(down_genes, file.path(output_dir, "Downregulated_genes_rebuilt.csv"), row.names = FALSE)

# 6. Variance-stabilized counts and normalized counts -----------
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)
write.csv(norm_counts, file.path(output_dir, "vst_normalized_counts_rebuilt.csv"))

# 7. Plots ----------------------------------------------------
png(file.path(output_dir, "MA_plot.png"), width = 1200, height = 900, res = 150)
plotMA(res, ylim = c(-5, 5), main = "Tumor vs Normal")
dev.off()

pca <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA: Tumor vs Normal")
ggsave(file.path(output_dir, "PCA_plot.png"), plot = pca, width = 7, height = 5, dpi = 300)

volcano_df <- res_df
volcano_df$significance <- "Not significant"
volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange > 1
] <- "Upregulated"
volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange < -1
] <- "Downregulated"

volcano_plot <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = -log10(pmax(padj, .Machine$double.xmin)),
      color = significance)
) +
  geom_point(alpha = 0.65, size = 1.2) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not significant" = "grey70"
  )) +
  labs(
    title = "Volcano Plot: Tumor vs Normal",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Result"
  ) +
  theme_classic()

ggsave(file.path(output_dir, "Volcano_plot.png"), volcano_plot,
       width = 7, height = 5, dpi = 300)

if (nrow(sig) > 1) {
  top_genes <- head(sig$ensembl_id[order(sig$padj)], min(30, nrow(sig)))
  annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
  pheatmap(
    assay(vsd)[top_genes, , drop = FALSE],
    scale = "row",
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    filename = file.path(output_dir, "Top_DEG_heatmap.png"),
    width = 8,
    height = 8
  )
}

# 8. Ensembl ID annotation ------------------------------------
# Remove version suffix (if present) before mapping.
ensembl_clean <- sub("\\..*$", "", sig$ensembl_id)
sig$ensembl_clean <- ensembl_clean

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensembl_clean),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
entrez_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensembl_clean),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

sig$gene_symbol <- unname(symbol_map[ensembl_clean])
sig$entrez_id <- unname(entrez_map[ensembl_clean])
write.csv(sig, file.path(output_dir, "DESeq2_significant_genes_annotated_rebuilt.csv"),
          row.names = FALSE)

# 9. GO Biological Process enrichment -------------------------
# Use mapped Entrez IDs; enrichment is run only if enough IDs map.
entrez_sig <- unique(na.omit(sig$entrez_id))
if (length(entrez_sig) >= 10) {
  ego <- enrichGO(
    gene = entrez_sig,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  go_df <- as.data.frame(ego)
  write.csv(go_df, file.path(output_dir, "GO_Biological_Process_Results_rebuilt.csv"),
            row.names = FALSE)
  if (nrow(go_df) > 0) {
    go_plot <- dotplot(ego, showCategory = min(15, nrow(go_df))) +
      ggtitle("GO Biological Process Enrichment")
    ggsave(file.path(output_dir, "GO_BP_dotplot.png"), go_plot,
           width = 8, height = 6, dpi = 300)
  }
} else {
  warning("Fewer than 10 significant genes mapped to Entrez IDs; GO enrichment skipped.")
}

# 10. KEGG enrichment -----------------------------------------
if (length(entrez_sig) >= 10) {
  ekegg <- enrichKEGG(
    gene = entrez_sig,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  kegg_df <- as.data.frame(ekegg)
  write.csv(kegg_df, file.path(output_dir, "KEGG_Pathway_Results_rebuilt.csv"),
            row.names = FALSE)
  if (nrow(kegg_df) > 0) {
    kegg_plot <- dotplot(ekegg, showCategory = min(15, nrow(kegg_df))) +
      ggtitle("KEGG Pathway Enrichment")
    ggsave(file.path(output_dir, "KEGG_dotplot.png"), kegg_plot,
           width = 8, height = 6, dpi = 300)
  }
} else {
  warning("Fewer than 10 significant genes mapped to Entrez IDs; KEGG enrichment skipped.")
}

# 11. Save session information and summary --------------------
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))

summary_lines <- c(
  "RNA-seq analysis: GSE81089 Tumor vs Normal",
  paste("Input file:", input_file),
  paste("Samples in raw count matrix:", ncol(count_matrix)),
  paste("Paired patients:", length(paired_patients)),
  paste("Paired samples:", ncol(count_matrix_paired)),
  paste("Genes after low-count filtering:", nrow(count_matrix_filtered)),
  paste("Significant DEGs (padj < 0.05 and |log2FC| > 1):", nrow(sig)),
  paste("Upregulated:", nrow(up_genes)),
  paste("Downregulated:", nrow(down_genes)),
  "",
  "Review sample metadata and design assumptions before reporting results."
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nFinished. Outputs saved in:\n", output_dir, "\n")
