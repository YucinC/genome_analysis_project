# ============================================================
# RNA-seq DESeq2 analysis using HTSeq-count output
# Reference: Canu assembly
# Comparison: Serum vs BH
#
# Improved version:
# 1. Merge DESeq2 results with Prokka annotation
# 2. MA plot with threshold lines, colors, and gene labels
# 3. PCA plot with stretched y-axis and sample labels
# 4. Sample distance heatmap
# 5. Volcano plot with threshold lines, colors, and gene labels
# 6. Summary tables for plots and labeled genes
# ============================================================


# -----------------------------
# 0. Load packages
# -----------------------------

# If packages are not installed, run these once:
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages(c("ggplot2", "ggrepel", "pheatmap", "dplyr", "readr"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(readr)
})


# -----------------------------
# 1. Set working paths
# -----------------------------

# Count files are searched under this folder.
# If your R working directory already contains RNA-Seq_BH and RNA-Seq_Serum folders,
# keep COUNT_DIR as ".".
COUNT_DIR <- "."

# Prokka annotation table.
# Use "/" instead of "\" in Windows paths.
ANNOTATION_FILE <- "A:/genome_project/E745_canu.tsv"

OUT_DIR <- file.path(dirname(COUNT_DIR), "deseq2_rna_canu_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Count directory:\n", COUNT_DIR, "\n")
cat("Annotation file:\n", ANNOTATION_FILE, "\n")
cat("Output directory:\n", OUT_DIR, "\n\n")

if (!file.exists(ANNOTATION_FILE)) {
  stop("Annotation file not found. Please check ANNOTATION_FILE.")
}


# -----------------------------
# 2. Find HTSeq-count files
# -----------------------------

count_files <- list.files(
  COUNT_DIR,
  pattern = "\\.counts\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

count_files <- sort(count_files)

if (length(count_files) == 0) {
  stop("No HTSeq-count files found. Please check COUNT_DIR.")
}

cat("Found count files:\n")
print(count_files)


# -----------------------------
# 3. Define sample information
# -----------------------------

sample_names <- basename(count_files)
sample_names <- sub("\\.canu\\.htseq\\.counts\\.txt$", "", sample_names)
sample_names <- sub("\\.htseq\\.counts\\.txt$", "", sample_names)
sample_names <- sub("\\.counts\\.txt$", "", sample_names)

# The dataset folder name should be RNA-Seq_BH or RNA-Seq_Serum.
dataset <- basename(dirname(count_files))

condition <- ifelse(
  dataset == "RNA-Seq_BH", "BH",
  ifelse(dataset == "RNA-Seq_Serum", "Serum", NA)
)

if (any(is.na(condition))) {
  stop("Some samples could not be assigned to BH or Serum. Check folder names.")
}

metadata <- data.frame(
  sample = sample_names,
  dataset = dataset,
  condition = condition,
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$sample

cat("\nSample metadata:\n")
print(metadata)

write.csv(
  metadata,
  file.path(OUT_DIR, "sample_metadata.csv"),
  quote = FALSE,
  row.names = TRUE
)


# -----------------------------
# 4. Read HTSeq-count files
# -----------------------------

read_htseq_count <- function(file) {
  df <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("gene_id", "count")
  )
  
  # Remove HTSeq summary rows:
  # __no_feature, __ambiguous, __too_low_aQual, __not_aligned, etc.
  df <- df[!grepl("^__", df$gene_id), ]
  
  df$count <- as.integer(df$count)
  return(df)
}

count_list <- lapply(count_files, read_htseq_count)
names(count_list) <- sample_names

# Check whether all files contain the same genes.
gene_sets <- lapply(count_list, function(x) x$gene_id)
all_genes <- Reduce(union, gene_sets)

cat("\nNumber of genes detected across all count files:", length(all_genes), "\n")

# Build count matrix.
count_matrix <- matrix(
  0L,
  nrow = length(all_genes),
  ncol = length(count_list),
  dimnames = list(all_genes, sample_names)
)

for (sample in names(count_list)) {
  df <- count_list[[sample]]
  count_matrix[df$gene_id, sample] <- df$count
}

# Reorder columns to match metadata.
count_matrix <- count_matrix[, rownames(metadata)]

write.csv(
  as.data.frame(count_matrix),
  file.path(OUT_DIR, "raw_count_matrix.csv"),
  quote = FALSE
)

cat("\nRaw count matrix dimension:\n")
print(dim(count_matrix))

cat("\nLibrary sizes based on CDS-assigned reads:\n")
print(colSums(count_matrix))


# -----------------------------
# 5. Prefilter low-count genes
# -----------------------------

# Keep genes with at least 10 counts across all samples.
keep <- rowSums(count_matrix) >= 10
count_matrix_filtered <- count_matrix[keep, ]

cat("\nGenes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")

write.csv(
  as.data.frame(count_matrix_filtered),
  file.path(OUT_DIR, "raw_count_matrix_filtered.csv"),
  quote = FALSE
)


# -----------------------------
# 6. Run DESeq2
# -----------------------------

metadata$condition <- factor(metadata$condition, levels = c("BH", "Serum"))

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix_filtered),
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

saveRDS(dds, file.path(OUT_DIR, "dds_rna_canu.rds"))

# Check size factors.
size_factor_table <- data.frame(
  sample = colnames(count_matrix_filtered),
  condition = metadata$condition,
  raw_assigned_counts = colSums(count_matrix_filtered),
  size_factor = sizeFactors(dds)
)

cat("\nDESeq2 size factors:\n")
print(size_factor_table)

write.csv(
  size_factor_table,
  file.path(OUT_DIR, "deseq2_size_factors.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Normalized counts.
normalized_counts <- counts(dds, normalized = TRUE)

write.csv(
  as.data.frame(normalized_counts),
  file.path(OUT_DIR, "normalized_counts.csv"),
  quote = FALSE
)


# -----------------------------
# 7. Differential expression result
# -----------------------------

# Serum vs BH:
# log2FoldChange > 0 means higher expression in Serum.
# log2FoldChange < 0 means higher expression in BH.
res <- results(dds, contrast = c("condition", "Serum", "BH"))
res <- res[order(res$padj), ]

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

write.csv(
  res_df,
  file.path(OUT_DIR, "Serum_vs_BH_all_results_unannotated.csv"),
  row.names = FALSE,
  quote = FALSE
)


# -----------------------------
# 8. Read and process Prokka annotation
# -----------------------------

annotation_raw <- read.delim(
  ANNOTATION_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("\nAnnotation table columns:\n")
print(colnames(annotation_raw))

# Prokka .tsv usually contains:
# locus_tag, ftype, length_bp, gene, EC_number, COG, product
#
# The matching key should usually be locus_tag if HTSeq was counted using Prokka GFF.
# This block is written to be tolerant of small column-name differences.

if (!("locus_tag" %in% colnames(annotation_raw))) {
  stop("The annotation file does not contain a 'locus_tag' column. Please check the Prokka TSV format.")
}

if (!("gene" %in% colnames(annotation_raw))) {
  annotation_raw$gene <- NA
}

if (!("product" %in% colnames(annotation_raw))) {
  annotation_raw$product <- NA
}

if (!("ftype" %in% colnames(annotation_raw))) {
  annotation_raw$ftype <- NA
}

if (!("COG" %in% colnames(annotation_raw))) {
  annotation_raw$COG <- NA
}

if (!("EC_number" %in% colnames(annotation_raw))) {
  annotation_raw$EC_number <- NA
}

annotation_df <- annotation_raw %>%
  mutate(
    gene_name = ifelse(
      !is.na(gene) & gene != "" & gene != "-",
      gene,
      NA
    ),
    product_name = ifelse(
      !is.na(product) & product != "" & product != "-",
      product,
      NA
    ),
    product_short = ifelse(
      !is.na(product_name),
      substr(product_name, 1, 45),
      NA
    ),
    plot_label = case_when(
      !is.na(gene_name) & !is.na(product_short) ~ paste0(gene_name, " / ", product_short),
      !is.na(gene_name) ~ gene_name,
      !is.na(product_short) ~ product_short,
      TRUE ~ locus_tag
    )
  ) %>%
  select(
    gene_id = locus_tag,
    ftype,
    gene_name,
    product_name,
    product_short,
    plot_label,
    COG,
    EC_number
  )

# If the same locus_tag appears more than once, keep the first.
annotation_df <- annotation_df[!duplicated(annotation_df$gene_id), ]

write.csv(
  annotation_df,
  file.path(OUT_DIR, "prokka_annotation_simplified.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Merge DESeq2 results with annotation.
res_annotated <- res_df %>%
  left_join(annotation_df, by = "gene_id") %>%
  mutate(
    gene_name = ifelse(is.na(gene_name), "", gene_name),
    product_name = ifelse(is.na(product_name), "", product_name),
    product_short = ifelse(is.na(product_short), "", product_short),
    plot_label = ifelse(is.na(plot_label) | plot_label == "", gene_id, plot_label),
    expression_change = case_when(
      !is.na(padj) & padj < 0.05 & !is.na(log2FoldChange) & log2FoldChange >= 1 ~ "Up in Serum",
      !is.na(padj) & padj < 0.05 & !is.na(log2FoldChange) & log2FoldChange <= -1 ~ "Higher in BH",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  res_annotated,
  file.path(OUT_DIR, "Serum_vs_BH_all_results_annotated.csv"),
  row.names = FALSE,
  quote = FALSE
)


# -----------------------------
# 9. Significant genes
# -----------------------------

sig <- res_annotated %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    !is.na(log2FoldChange),
    abs(log2FoldChange) >= 1
  )

up <- sig %>%
  filter(log2FoldChange > 0)

down <- sig %>%
  filter(log2FoldChange < 0)

write.csv(
  sig,
  file.path(OUT_DIR, "Serum_vs_BH_significant_padj0.05_log2FC1_annotated.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  up,
  file.path(OUT_DIR, "Serum_vs_BH_up_in_serum_annotated.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  down,
  file.path(OUT_DIR, "Serum_vs_BH_higher_in_BH_annotated.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nDESeq2 result summary:\n")
print(summary(res))

cat("\nNumber of significant genes, padj < 0.05 and |log2FC| >= 1:", nrow(sig), "\n")
cat("Upregulated in Serum:", nrow(up), "\n")
cat("Higher expression in BH:", nrow(down), "\n")


# -----------------------------
# 10. Select genes to label in MA and volcano plots
# -----------------------------

# Label strategy:
# - Top 5 significant genes upregulated in Serum by adjusted p-value
# - Top 5 significant genes higher in BH by adjusted p-value
# This keeps the figure readable and biologically directional.

top_up_for_label <- up %>%
  arrange(padj) %>%
  slice_head(n = 5)

top_down_for_label <- down %>%
  arrange(padj) %>%
  slice_head(n = 5)

label_genes <- bind_rows(top_up_for_label, top_down_for_label) %>%
  distinct(gene_id, .keep_all = TRUE)

write.csv(
  label_genes,
  file.path(OUT_DIR, "Serum_vs_BH_labeled_genes_in_plots.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nGenes selected for plot labels:\n")
print(label_genes[, c("gene_id", "gene_name", "product_short", "log2FoldChange", "padj", "expression_change")])


# -----------------------------
# 11. VST transformation
# -----------------------------

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(OUT_DIR, "vst_rna_canu.rds"))


# -----------------------------
# 12. PCA plot with sample labels
# -----------------------------

pca_matrix <- t(assay(vsd))
pca <- prcomp(pca_matrix)

percent_var <- pca$sdev^2 / sum(pca$sdev^2)
percent_var <- round(percent_var * 100, 2)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = metadata[rownames(pca$x), "condition"]
)

write.csv(
  pca_df,
  file.path(OUT_DIR, "PCA_coordinates.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Stretch the y-axis visually by using a larger figure height and explicit expansion.
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    title = "PCA plot of RNA-seq samples",
    subtitle = "VST-transformed counts; sample labels are shown",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))

ggsave(
  filename = file.path(OUT_DIR, "PCA_plot_Serum_vs_BH_labeled.pdf"),
  plot = p_pca,
  width = 7,
  height = 7
)

ggsave(
  filename = file.path(OUT_DIR, "PCA_plot_Serum_vs_BH_labeled.png"),
  plot = p_pca,
  width = 7,
  height = 7,
  dpi = 300
)


# -----------------------------
# 13. Sample distance heatmap
# -----------------------------

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- colnames(vsd)

annotation_col <- data.frame(
  condition = metadata[colnames(vsd), "condition"]
)
rownames(annotation_col) <- colnames(vsd)

pdf(
  file.path(OUT_DIR, "Sample_distance_heatmap.pdf"),
  width = 7,
  height = 6
)

pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = annotation_col,
  main = "Sample-to-sample distance based on VST counts",
  fontsize = 11,
  fontsize_row = 10,
  fontsize_col = 10
)

dev.off()

png(
  file.path(OUT_DIR, "Sample_distance_heatmap.png"),
  width = 2100,
  height = 1800,
  res = 300
)

pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = annotation_col,
  main = "Sample-to-sample distance based on VST counts",
  fontsize = 11,
  fontsize_row = 10,
  fontsize_col = 10
)

dev.off()

write.csv(
  sample_dist_matrix,
  file.path(OUT_DIR, "Sample_distance_matrix.csv"),
  quote = FALSE
)


# -----------------------------
# 14. MA plot with clearer annotation labels
# -----------------------------

ma_df <- res_annotated %>%
  mutate(
    baseMean_for_plot = ifelse(baseMean <= 0 | is.na(baseMean), NA, baseMean),
    label_this_gene = gene_id %in% label_genes$gene_id,
    expression_change = factor(
      expression_change,
      levels = c("Higher in BH", "Not significant", "Up in Serum")
    )
  ) %>%
  filter(!is.na(baseMean_for_plot), !is.na(log2FoldChange))

p_ma <- ggplot(ma_df, aes(x = baseMean_for_plot, y = log2FoldChange, color = expression_change)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5, color = "grey40") +
  
  # Use label boxes instead of plain text so labels stand out
  geom_label_repel(
    data = ma_df %>% filter(label_this_gene),
    aes(label = plot_label),
    color = "black",             # label text color
    fill = "white",              # label background
    label.size = 0.25,           # border thickness of label box
    size = 3.0,
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey35",
    show.legend = FALSE
  ) +
  
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Higher in BH" = "#2C7FB8",
      "Not significant" = "grey70",
      "Up in Serum" = "#D95F0E"
    )
  ) +
  labs(
    title = "MA plot: Serum vs BH",
    subtitle = paste0(
      "Significant genes: padj < 0.05 and |log2FC| >= 1\n",
      "Up in Serum = ", nrow(up), "; Higher in BH = ", nrow(down)
    ),
    x = "Mean of normalized counts",
    y = "log2 fold change: Serum / BH",
    color = "Expression change"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

ggsave(
  filename = file.path(OUT_DIR, "MA_plot_Serum_vs_BH_annotated.pdf"),
  plot = p_ma,
  width = 8.5,
  height = 6
)

ggsave(
  filename = file.path(OUT_DIR, "MA_plot_Serum_vs_BH_annotated.png"),
  plot = p_ma,
  width = 8.5,
  height = 6,
  dpi = 300
)


# -----------------------------
# 15. Volcano plot with clearer annotation labels
# -----------------------------

volcano_df <- res_annotated %>%
  mutate(
    padj_for_plot = case_when(
      is.na(padj) ~ NA_real_,
      padj == 0 ~ .Machine$double.xmin,
      TRUE ~ padj
    ),
    neg_log10_padj = -log10(padj_for_plot),
    label_this_gene = gene_id %in% label_genes$gene_id,
    expression_change = factor(
      expression_change,
      levels = c("Higher in BH", "Not significant", "Up in Serum")
    )
  ) %>%
  filter(!is.na(log2FoldChange), !is.na(neg_log10_padj))

p_volcano <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = neg_log10_padj, color = expression_change)
) +
  geom_point(alpha = 0.75, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5, color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "grey40") +
  
  # Clear label boxes
  geom_label_repel(
    data = volcano_df %>% filter(label_this_gene),
    aes(label = plot_label),
    color = "black",             # label text color
    fill = "white",              # label background
    label.size = 0.25,
    size = 3.0,
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey35",
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Higher in BH" = "#2C7FB8",
      "Not significant" = "grey70",
      "Up in Serum" = "#D95F0E"
    )
  ) +
  labs(
    title = "Volcano plot: Serum vs BH",
    subtitle = paste0(
      "Right: upregulated in Serum; Left: higher expression in BH\n",
      "Threshold: padj < 0.05 and |log2FC| >= 1"
    ),
    x = "log2 fold change: Serum / BH",
    y = "-log10 adjusted p-value",
    color = "Expression change"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

ggsave(
  filename = file.path(OUT_DIR, "Volcano_plot_Serum_vs_BH_annotated.pdf"),
  plot = p_volcano,
  width = 8.5,
  height = 6
)

ggsave(
  filename = file.path(OUT_DIR, "Volcano_plot_Serum_vs_BH_annotated.png"),
  plot = p_volcano,
  width = 8.5,
  height = 6,
  dpi = 300
)


# -----------------------------
# 16. Save top genes
# -----------------------------

top_20 <- res_annotated %>%
  arrange(padj) %>%
  slice_head(n = 20)

write.csv(
  top_20,
  file.path(OUT_DIR, "Serum_vs_BH_top20_by_padj_annotated.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nTop 20 genes by adjusted p-value:\n")
print(top_20[, c("gene_id", "gene_name", "product_short", "baseMean", "log2FoldChange", "padj", "expression_change")])


# -----------------------------
# 17. Plot summary table
# -----------------------------

plot_summary <- data.frame(
  output_file = c(
    "MA_plot_Serum_vs_BH_annotated.pdf/png",
    "PCA_plot_Serum_vs_BH_labeled.pdf/png",
    "Sample_distance_heatmap.pdf/png",
    "Volcano_plot_Serum_vs_BH_annotated.pdf/png",
    "Serum_vs_BH_all_results_annotated.csv",
    "Serum_vs_BH_labeled_genes_in_plots.csv"
  ),
  result_type = c(
    "MA plot",
    "PCA plot",
    "Sample distance heatmap",
    "Volcano plot",
    "Annotated DESeq2 result table",
    "Labeled gene table"
  ),
  main_question = c(
    "Which genes show expression changes across different average expression levels?",
    "Do RNA-seq samples separate by condition, and are replicates consistent?",
    "Which samples are globally similar or different based on VST-transformed expression?",
    "Which genes show both large fold change and strong statistical significance?",
    "What are the differential expression statistics and gene annotations for all genes?",
    "Which genes were labeled in MA and volcano plots?"
  ),
  how_to_read = c(
    "x-axis is mean normalized expression; y-axis is log2FC. Positive log2FC means higher in Serum, negative means higher in BH.",
    "Each point is one sample. Separation along PC1 or PC2 indicates major expression differences between samples.",
    "Rows and columns are samples. Smaller distance means more similar global expression profiles.",
    "x-axis is log2FC; y-axis is -log10 adjusted p-value. Upper-right genes are significantly upregulated in Serum; upper-left genes are higher in BH.",
    "Use baseMean, log2FoldChange, padj, gene_name, and product_name to interpret each gene.",
    "These are the top significant genes selected for direct labels in the plots."
  ),
  what_to_check = c(
    "Check whether significant genes are concentrated at very low expression, whether there is global bias, and whether extreme outliers exist.",
    "Check whether BH and Serum samples separate clearly and whether samples from the same condition cluster together.",
    "Check whether biological replicates are closer to each other than to the other condition, and whether any sample is an outlier.",
    "Check the balance of upregulated and downregulated genes, the top labeled genes, and whether extreme fold changes are biologically plausible.",
    "Check top genes by padj and whether their annotated functions match the expected serum-response biology.",
    "Check whether the labels are readable and whether the selected genes are suitable examples for interpretation."
  ),
  current_result_summary = c(
    paste0("Significant genes = ", nrow(sig), "; Up in Serum = ", nrow(up), "; Higher in BH = ", nrow(down), "."),
    paste0("PC1 explains ", percent_var[1], "% variance; PC2 explains ", percent_var[2], "% variance."),
    "Use this heatmap as a supplement to PCA, especially when PC2 explains very little variance.",
    paste0("Significant genes = ", nrow(sig), "; Up in Serum = ", nrow(up), "; Higher in BH = ", nrow(down), "."),
    paste0("Annotated result table contains ", nrow(res_annotated), " genes after DESeq2 filtering."),
    paste0("Number of labeled genes = ", nrow(label_genes), ".")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  plot_summary,
  file.path(OUT_DIR, "Serum_vs_BH_plot_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  c(
    "RNA-seq plot interpretation summary",
    "===================================",
    "",
    paste0("Comparison: Serum vs BH"),
    paste0("Positive log2FoldChange: higher expression in Serum"),
    paste0("Negative log2FoldChange: higher expression in BH"),
    "",
    paste0("Total genes tested after filtering: ", nrow(res_annotated)),
    paste0("Significant genes, padj < 0.05 and |log2FC| >= 1: ", nrow(sig)),
    paste0("Upregulated in Serum: ", nrow(up)),
    paste0("Higher expression in BH: ", nrow(down)),
    "",
    "MA plot:",
    "The MA plot shows whether differential expression depends on the average expression level.",
    "Genes above zero have higher expression in Serum, while genes below zero have higher expression in BH.",
    "The dashed horizontal lines mark log2FC = +1 and -1.",
    "",
    "PCA plot:",
    "The PCA plot shows whether samples separate by condition based on global expression profiles.",
    paste0("PC1 explains ", percent_var[1], "% of the variance."),
    paste0("PC2 explains ", percent_var[2], "% of the variance."),
    "If PC2 is very small, the sample distance heatmap should be used as a complementary sample-level QC plot.",
    "",
    "Sample distance heatmap:",
    "This heatmap shows pairwise distances between samples after VST transformation.",
    "Replicates from the same condition should generally be closer to each other.",
    "",
    "Volcano plot:",
    "The volcano plot highlights genes with both large fold changes and strong statistical significance.",
    "The right side shows genes upregulated in Serum, while the left side shows genes with higher expression in BH.",
    "The vertical dashed lines mark log2FC = +1 and -1, and the horizontal dashed line marks padj = 0.05."
  ),
  con = file.path(OUT_DIR, "Serum_vs_BH_plot_interpretation.txt")
)


# -----------------------------
# 18. Finish
# -----------------------------

cat("\nAnalysis finished successfully.\n")
cat("Results saved to:\n", OUT_DIR, "\n\n")

cat("Main improved output files:\n")
cat("- Serum_vs_BH_all_results_annotated.csv\n")
cat("- Serum_vs_BH_significant_padj0.05_log2FC1_annotated.csv\n")
cat("- Serum_vs_BH_up_in_serum_annotated.csv\n")
cat("- Serum_vs_BH_higher_in_BH_annotated.csv\n")
cat("- Serum_vs_BH_labeled_genes_in_plots.csv\n")
cat("- Serum_vs_BH_plot_summary.csv\n")
cat("- Serum_vs_BH_plot_interpretation.txt\n")
cat("- PCA_plot_Serum_vs_BH_labeled.pdf/png\n")
cat("- Sample_distance_heatmap.pdf/png\n")
cat("- MA_plot_Serum_vs_BH_annotated.pdf/png\n")
cat("- Volcano_plot_Serum_vs_BH_annotated.pdf/png\n")

