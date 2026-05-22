# ============================================================
# Tn-seq window-to-gene aggregated DESeq2 analysis
# Reference: Canu assembly
#
# Input:
#   1. tnseq_gene_counts_from_25bp_windows.matrix.tsv
#   2. tnseq_sample_metadata.tsv
#
# Meaning:
#   25 bp window counts were overlapped with CDS regions and
#   aggregated to gene-level counts.
#
# Output:
#   Annotated DESeq2 tables, PCA, sample-distance heatmap,
#   annotated MA plots, annotated volcano plots, and result summaries.
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
})

# -----------------------------
# 1. Basic settings
# -----------------------------

COUNT_MATRIX_FILE <- "../tnseq_counting/tnseq_gene_counts_from_25bp_windows.matrix.tsv"
METADATA_FILE <- "../tnseq_counting/tnseq_sample_metadata.tsv"


COUNT_MATRIX_FILE <- "../tnseq_counting/tnseq_gene_counts_from_25bp_windows.matrix.tsv"
METADATA_FILE <- "../tnseq_counting/tnseq_sample_metadata.tsv"
ANNOTATION_FILE <- "../E745_canu.tsv"

OUT_DIR <- "deseq_window_gene_results"
PLOT_DIR <- file.path(OUT_DIR, "plots")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

conditions <- c("Tn-Seq_BHI", "Tn-Seq_HSerum", "Tn-Seq_Serum")

LFC_CUTOFF <- 1
PADJ_CUTOFF <- 0.05
N_LABEL_EACH_DIRECTION <- 5

cat("Count matrix file:", COUNT_MATRIX_FILE, "\n")
cat("Metadata file:", METADATA_FILE, "\n")
cat("Annotation file:", ANNOTATION_FILE, "\n")
cat("Output directory:", OUT_DIR, "\n\n")


# -----------------------------
# 2. Helper functions
# -----------------------------

clean_feature_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^ID=", "", x)
  x <- sub("^gene-", "", x)
  x <- sub("^cds-", "", x)
  x <- sub("^CDS:", "", x)
  x <- sub(";.*$", "", x)
  return(x)
}

clean_sample_name <- function(x) {
  x <- basename(x)
  x <- sub("\\.25bp_window\\.counts\\.tsv$", "", x)
  x <- sub("\\.canu\\.htseq\\.counts\\.txt$", "", x)
  x <- sub("\\.htseq\\.counts\\.txt$", "", x)
  x <- sub("\\.counts\\.txt$", "", x)
  x <- sub("\\.txt$", "", x)
  x <- sub("_pass$", "", x)
  return(x)
}

safe_pdf_save <- function(filename, plot, width, height) {
  try(
    ggsave(filename = filename, plot = plot, width = width, height = height),
    silent = TRUE
  )
}

read_prokka_annotation <- function(annotation_file) {
  if (!file.exists(annotation_file)) {
    warning("Annotation file not found. Plot labels will fall back to feature IDs.")
    return(NULL)
  }
  
  anno <- read.delim(
    annotation_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
  
  message("Annotation columns detected:")
  message(paste(colnames(anno), collapse = ", "))
  
  if (!("locus_tag" %in% colnames(anno))) {
    stop("The Prokka annotation file does not contain a 'locus_tag' column.")
  }
  
  if (!("gene" %in% colnames(anno))) anno$gene <- NA
  if (!("product" %in% colnames(anno))) anno$product <- NA
  if (!("ftype" %in% colnames(anno))) anno$ftype <- NA
  if (!("COG" %in% colnames(anno))) anno$COG <- NA
  if (!("EC_number" %in% colnames(anno))) anno$EC_number <- NA
  
  anno2 <- anno %>%
    mutate(
      feature = clean_feature_id(locus_tag),
      
      gene_name = ifelse(
        is.na(gene) | gene == "" | gene == "-",
        NA,
        gene
      ),
      
      product_name = ifelse(
        is.na(product) | product == "" | product == "-",
        NA,
        product
      ),
      
      product_short = ifelse(
        !is.na(product_name),
        substr(product_name, 1, 55),
        NA
      ),
      
      has_gene_name = !is.na(gene_name),
      has_product_name = !is.na(product_name),
      
      annotation_priority = case_when(
        ftype == "CDS" & has_gene_name & has_product_name ~ 1,
        ftype == "CDS" & has_gene_name ~ 2,
        ftype == "CDS" & has_product_name ~ 3,
        ftype == "CDS" ~ 4,
        has_gene_name & has_product_name ~ 5,
        has_gene_name ~ 6,
        has_product_name ~ 7,
        TRUE ~ 8
      )
    ) %>%
    arrange(feature, annotation_priority) %>%
    group_by(feature) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      gene_only_label = case_when(
        !is.na(gene_name) & gene_name != "" ~ gene_name,
        TRUE ~ feature
      ),
      
      gene_product_label = case_when(
        !is.na(gene_name) & gene_name != "" &
          !is.na(product_short) & product_short != "" ~ paste0(gene_name, " / ", product_short),
        
        !is.na(gene_name) & gene_name != "" ~ gene_name,
        
        !is.na(product_short) & product_short != "" ~ product_short,
        
        TRUE ~ feature
      )
    ) %>%
    select(
      feature,
      ftype,
      gene_name,
      product_name,
      product_short,
      gene_only_label,
      gene_product_label,
      COG,
      EC_number
    )
  
  return(anno2)
}

annotate_result <- function(res_df, annotation_df) {
  res_df$feature <- clean_feature_id(res_df$feature)

  if (is.null(annotation_df)) {
    res_df$ftype <- NA
    res_df$gene_name <- NA
    res_df$product_name <- NA
    res_df$product_short <- NA
    res_df$plot_label <- res_df$feature
    res_df$COG <- NA
    res_df$EC_number <- NA
    res_df$annotation_status <- "no_annotation_file"
    return(res_df)
  }

  annotation_df$feature <- clean_feature_id(annotation_df$feature)

  res_annotated <- res_df %>%
    left_join(annotation_df, by = "feature") %>%
    mutate(
      annotation_status = ifelse(
        is.na(product_name) & is.na(gene_name),
        "not_annotated_or_no_product",
        "annotated"
      ),
      plot_label = case_when(
        !is.na(gene_name) & gene_name != "" &
          !is.na(product_short) & product_short != "" ~ paste0(gene_name, " / ", product_short),
        !is.na(gene_name) & gene_name != "" ~ gene_name,
        !is.na(product_short) & product_short != "" ~ product_short,
        !is.na(product_name) & product_name != "" ~ substr(product_name, 1, 55),
        TRUE ~ feature
      )
    )

  n_total <- nrow(res_annotated)
  n_annotated <- sum(res_annotated$annotation_status == "annotated", na.rm = TRUE)
  message("Annotation matching summary:")
  message("  Total result features: ", n_total)
  message("  Annotated with gene/product: ", n_annotated)
  message("  Annotation rate: ", round(100 * n_annotated / n_total, 2), "%")

  return(res_annotated)
}

classify_tnseq_change <- function(df, lfc_cutoff = 1, padj_cutoff = 0.05) {
  df$tnseq_change <- "Not significant"

  higher_idx <- !is.na(df$padj) &
    df$padj < padj_cutoff &
    !is.na(df$log2FoldChange) &
    df$log2FoldChange >= lfc_cutoff

  lower_idx <- !is.na(df$padj) &
    df$padj < padj_cutoff &
    !is.na(df$log2FoldChange) &
    df$log2FoldChange <= -lfc_cutoff

  df$tnseq_change[higher_idx] <- "Higher mutant abundance"
  df$tnseq_change[lower_idx] <- "Lower mutant abundance"

  df$tnseq_change <- factor(
    df$tnseq_change,
    levels = c("Lower mutant abundance", "Not significant", "Higher mutant abundance")
  )

  return(df)
}

select_label_features <- function(df, n_each_direction = 5, lfc_cutoff = 1, padj_cutoff = 0.05) {
  sig_df <- df %>%
    filter(
      !is.na(padj),
      padj < padj_cutoff,
      !is.na(log2FoldChange),
      abs(log2FoldChange) >= lfc_cutoff
    )

  top_higher <- sig_df %>%
    filter(log2FoldChange > 0) %>%
    arrange(padj) %>%
    slice_head(n = n_each_direction)

  top_lower <- sig_df %>%
    filter(log2FoldChange < 0) %>%
    arrange(padj) %>%
    slice_head(n = n_each_direction)

  bind_rows(top_higher, top_lower) %>%
    distinct(feature, .keep_all = TRUE)
}

comparison_label <- function(comparison_name) {
  x <- gsub("_vs_", " / ", comparison_name)
  x <- gsub("Tn-Seq_", "", x)
  return(x)
}

detect_metadata_columns <- function(metadata) {
  sample_candidates <- c("sample", "Sample", "sample_id", "SampleID", "id", "ID")
  condition_candidates <- c("condition", "Condition", "dataset", "Dataset", "group", "Group")

  sample_col <- intersect(sample_candidates, colnames(metadata))[1]
  condition_col <- intersect(condition_candidates, colnames(metadata))[1]

  if (is.na(sample_col)) {
    stop("Could not detect sample column in metadata file.")
  }

  if (is.na(condition_col)) {
    stop("Could not detect condition column in metadata file.")
  }

  return(list(sample_col = sample_col, condition_col = condition_col))
}


# -----------------------------
# 3. Read count matrix and metadata
# -----------------------------

if (!file.exists(COUNT_MATRIX_FILE)) {
  stop("COUNT_MATRIX_FILE not found. Please check the file path.")
}

if (!file.exists(METADATA_FILE)) {
  stop("METADATA_FILE not found. Please check the file path.")
}

count_raw <- read.delim(
  COUNT_MATRIX_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(count_raw) < 2) {
  stop("Count matrix has fewer than 2 columns. Please check file format.")
}

feature_col <- colnames(count_raw)[1]
count_features <- clean_feature_id(count_raw[[feature_col]])

count_matrix <- count_raw[, -1, drop = FALSE]
rownames(count_matrix) <- count_features
colnames(count_matrix) <- clean_sample_name(colnames(count_matrix))

count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric"
count_matrix[is.na(count_matrix)] <- 0
count_matrix <- round(count_matrix)

metadata_raw <- read.delim(
  METADATA_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata_cols <- detect_metadata_columns(metadata_raw)

sample_info <- data.frame(
  sample = clean_sample_name(metadata_raw[[metadata_cols$sample_col]]),
  condition = metadata_raw[[metadata_cols$condition_col]],
  stringsAsFactors = FALSE
)

sample_info$condition <- factor(sample_info$condition, levels = conditions)

if (any(is.na(sample_info$condition))) {
  stop("Some metadata conditions are NA after factor conversion. Please check condition names.")
}

common_samples <- intersect(sample_info$sample, colnames(count_matrix))

if (length(common_samples) == 0) {
  stop("No common sample names between count matrix and metadata.")
}

sample_info <- sample_info %>% filter(sample %in% common_samples)
count_matrix <- count_matrix[, sample_info$sample, drop = FALSE]

write.csv(sample_info, file.path(OUT_DIR, "sample_metadata.csv"), row.names = FALSE, quote = FALSE)
write.csv(as.data.frame(count_matrix), file.path(OUT_DIR, "window_gene_count_matrix.csv"), quote = FALSE)

cat("Count matrix dimension:\n")
print(dim(count_matrix))
cat("\nSample information:\n")
print(sample_info)


# -----------------------------
# 4. DESeq2
# -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

features_before_filtering <- nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10, ]
features_after_filtering <- nrow(dds)

dds <- DESeq(dds)

saveRDS(dds, file.path(OUT_DIR, "dds_window_gene_tnseq.rds"))

norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts), file.path(OUT_DIR, "window_gene_normalized_counts.csv"), quote = FALSE)

size_factor_table <- data.frame(
  sample = colnames(dds),
  condition = colData(dds)$condition,
  raw_counts_after_filtering = colSums(counts(dds, normalized = FALSE)),
  size_factor = sizeFactors(dds),
  stringsAsFactors = FALSE
)

write.csv(size_factor_table, file.path(OUT_DIR, "deseq2_size_factors.csv"), row.names = FALSE, quote = FALSE)


# -----------------------------
# 5. Annotation
# -----------------------------

annotation_df <- read_prokka_annotation(ANNOTATION_FILE)

if (!is.null(annotation_df)) {
  write.csv(annotation_df, file.path(OUT_DIR, "prokka_annotation_simplified.csv"), row.names = FALSE, quote = FALSE)
}


# -----------------------------
# 6. Global QC plots
# -----------------------------

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(OUT_DIR, "vst_window_gene_tnseq.rds"))

pca <- prcomp(t(assay(vsd)))
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = colData(vsd)[rownames(pca$x), "condition"],
  stringsAsFactors = FALSE
)

write.csv(pca_df, file.path(OUT_DIR, "PCA_coordinates.csv"), row.names = FALSE, quote = FALSE)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3.5, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.35, 0.35))) +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.20))) +
  labs(
    title = "PCA plot of window-to-gene Tn-seq samples",
    subtitle = "VST-transformed gene-level counts aggregated from 25 bp windows",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        legend.position = "right")

ggsave(file.path(PLOT_DIR, "PCA_plot_window_gene_Tnseq_labeled.png"), p_pca, width = 7.5, height = 7, dpi = 300)
safe_pdf_save(file.path(PLOT_DIR, "PCA_plot_window_gene_Tnseq_labeled.pdf"), p_pca, width = 7.5, height = 7)

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

annotation_col <- data.frame(condition = colData(vsd)$condition)
rownames(annotation_col) <- colnames(vsd)

write.csv(sample_dist_matrix, file.path(OUT_DIR, "sample_distance_matrix.csv"), quote = FALSE)

png(file.path(PLOT_DIR, "Sample_distance_heatmap_window_gene_Tnseq.png"), width = 2100, height = 1800, res = 300)
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

try({
  pdf(file.path(PLOT_DIR, "Sample_distance_heatmap_window_gene_Tnseq.pdf"), width = 7, height = 6)
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
}, silent = TRUE)


# -----------------------------
# 7. Pairwise comparisons and plots
# -----------------------------

comparisons <- list(
  Serum_vs_BHI = c("condition", "Tn-Seq_Serum", "Tn-Seq_BHI"),
  HSerum_vs_BHI = c("condition", "Tn-Seq_HSerum", "Tn-Seq_BHI"),
  Serum_vs_HSerum = c("condition", "Tn-Seq_Serum", "Tn-Seq_HSerum")
)

summary_rows <- list()

for (name in names(comparisons)) {
  cat("\nRunning comparison:", name, "\n")

  res <- results(dds, contrast = comparisons[[name]])
  res <- res[order(res$padj), ]

  res_df <- as.data.frame(res)
  res_df$feature <- rownames(res_df)
  res_df <- res_df[, c("feature", setdiff(colnames(res_df), "feature"))]

  res_df <- annotate_result(res_df, annotation_df)
  res_df <- classify_tnseq_change(res_df, LFC_CUTOFF, PADJ_CUTOFF)

  annotation_check <- res_df %>% count(annotation_status)
  write.csv(annotation_check, file.path(OUT_DIR, paste0(name, "_annotation_matching_summary.csv")), row.names = FALSE, quote = FALSE)

  write.csv(
    res_df,
    file.path(OUT_DIR, paste0(name, "_window_gene_DESeq2_results_annotated.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  sig_padj <- res_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF)
  sig_lfc <- res_df %>% filter(!is.na(padj), padj < PADJ_CUTOFF, !is.na(log2FoldChange), abs(log2FoldChange) >= LFC_CUTOFF)

  write.csv(sig_padj, file.path(OUT_DIR, paste0(name, "_window_gene_significant_padj0.05.csv")), row.names = FALSE, quote = FALSE)
  write.csv(sig_lfc, file.path(OUT_DIR, paste0(name, "_window_gene_significant_padj0.05_log2FC1.csv")), row.names = FALSE, quote = FALSE)

  label_df <- select_label_features(res_df, N_LABEL_EACH_DIRECTION, LFC_CUTOFF, PADJ_CUTOFF)
  write.csv(label_df, file.path(OUT_DIR, paste0(name, "_window_gene_labeled_features.csv")), row.names = FALSE, quote = FALSE)

  n_higher <- sum(res_df$tnseq_change == "Higher mutant abundance", na.rm = TRUE)
  n_lower <- sum(res_df$tnseq_change == "Lower mutant abundance", na.rm = TRUE)

  ma_df <- res_df %>%
    mutate(
      baseMean_for_plot = ifelse(baseMean <= 0 | is.na(baseMean), NA, baseMean),
      label_this_feature = feature %in% label_df$feature
    ) %>%
    filter(!is.na(baseMean_for_plot), !is.na(log2FoldChange))

  p_ma <- ggplot(ma_df, aes(x = baseMean_for_plot, y = log2FoldChange, color = tnseq_change)) +
    geom_point(alpha = 0.7, size = 1.4) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_label_repel(
      data = ma_df %>% filter(label_this_feature),
      aes(label = plot_label),
      color = "black",
      fill = "white",
      label.size = 0.25,
      size = 3.0,
      max.overlaps = Inf,
      box.padding = 0.45,
      point.padding = 0.30,
      min.segment.length = 0,
      segment.color = "grey35",
      show.legend = FALSE
    ) +
    scale_x_log10() +
    scale_color_manual(values = c("Lower mutant abundance" = "#2C7FB8", "Not significant" = "grey70", "Higher mutant abundance" = "#D95F0E")) +
    labs(
      title = paste0("MA plot: ", name),
      subtitle = paste0("Significant features: padj < ", PADJ_CUTOFF, " and |log2FC| >= ", LFC_CUTOFF, "\n",
                        "Higher mutant abundance = ", n_higher, "; Lower mutant abundance = ", n_lower),
      x = "Mean of normalized insertion counts",
      y = paste0("log2 fold change: ", comparison_label(name)),
      color = "Tn-seq change"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 9),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.position = "right")

  ggsave(file.path(PLOT_DIR, paste0(name, "_window_gene_MA_plot_annotated.png")), p_ma, width = 8.5, height = 6, dpi = 300)
  safe_pdf_save(file.path(PLOT_DIR, paste0(name, "_window_gene_MA_plot_annotated.pdf")), p_ma, width = 8.5, height = 6)

  volcano_df <- res_df %>%
    mutate(
      padj_for_plot = case_when(is.na(padj) ~ NA_real_, padj == 0 ~ .Machine$double.xmin, TRUE ~ padj),
      neg_log10_padj = -log10(padj_for_plot),
      label_this_feature = feature %in% label_df$feature
    ) %>%
    filter(!is.na(log2FoldChange), !is.na(neg_log10_padj))

  p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = tnseq_change)) +
    geom_point(alpha = 0.75, size = 1.4) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_label_repel(
      data = volcano_df %>% filter(label_this_feature),
      aes(label = plot_label),
      color = "black",
      fill = "white",
      label.size = 0.25,
      size = 3.0,
      max.overlaps = Inf,
      box.padding = 0.45,
      point.padding = 0.30,
      min.segment.length = 0,
      segment.color = "grey35",
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("Lower mutant abundance" = "#2C7FB8", "Not significant" = "grey70", "Higher mutant abundance" = "#D95F0E")) +
    labs(
      title = paste0("Volcano plot: ", name),
      subtitle = paste0("Right: higher mutant abundance in numerator condition; left: lower mutant abundance\n",
                        "Threshold: padj < ", PADJ_CUTOFF, " and |log2FC| >= ", LFC_CUTOFF),
      x = paste0("log2 fold change: ", comparison_label(name)),
      y = "-log10 adjusted p-value",
      color = "Tn-seq change"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 9),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.position = "right")

  ggsave(file.path(PLOT_DIR, paste0(name, "_window_gene_volcano_plot_annotated.png")), p_volcano, width = 8.5, height = 6, dpi = 300)
  safe_pdf_save(file.path(PLOT_DIR, paste0(name, "_window_gene_volcano_plot_annotated.pdf")), p_volcano, width = 8.5, height = 6)

  top20 <- res_df %>% filter(!is.na(padj)) %>% arrange(padj) %>% slice_head(n = 20)
  write.csv(top20, file.path(OUT_DIR, paste0(name, "_window_gene_top20_by_padj_annotated.csv")), row.names = FALSE, quote = FALSE)

  summary_rows[[name]] <- data.frame(
    analysis = "window_to_gene_aggregated",
    comparison = name,
    tested_features = nrow(res_df),
    annotated_features = sum(res_df$annotation_status == "annotated", na.rm = TRUE),
    significant_padj_0.05 = nrow(sig_padj),
    significant_padj_0.05_log2FC_1 = nrow(sig_lfc),
    higher_mutant_abundance = n_higher,
    lower_mutant_abundance = n_lower,
    labeled_features = nrow(label_df),
    top_feature = ifelse(nrow(top20) > 0, top20$feature[1], NA),
    top_label = ifelse(nrow(top20) > 0, top20$plot_label[1], NA),
    top_log2FC = ifelse(nrow(top20) > 0, top20$log2FoldChange[1], NA),
    top_padj = ifelse(nrow(top20) > 0, top20$padj[1], NA),
    stringsAsFactors = FALSE
  )
}


# -----------------------------
# 8. Final summary
# -----------------------------

summary_table <- bind_rows(summary_rows)

write.csv(summary_table, file.path(OUT_DIR, "window_gene_Tnseq_DESeq2_result_summary.csv"), row.names = FALSE, quote = FALSE)

writeLines(
  c(
    "Tn-seq DESeq2 result summary: window-to-gene aggregated analysis",
    "=================================================================",
    "",
    "Interpretation:",
    "This analysis uses gene-level counts aggregated from 25 bp genomic windows overlapping CDS regions.",
    "It is not RNA-seq expression analysis.",
    "Positive log2FoldChange means higher mutant/insertion abundance in the numerator condition.",
    "Negative log2FoldChange means lower mutant/insertion abundance in the numerator condition.",
    "A significantly lower mutant abundance may suggest that disruption of the corresponding gene reduces fitness in the numerator condition.",
    "",
    paste0("Features before filtering: ", features_before_filtering),
    paste0("Features after filtering: ", features_after_filtering),
    "",
    "PCA and sample distance heatmap:",
    paste0("PC1 explains ", percent_var[1], "% variance."),
    paste0("PC2 explains ", percent_var[2], "% variance."),
    "Use these plots to check whether replicates cluster by condition and whether any sample is an outlier.",
    "",
    "Pairwise comparison summary:",
    capture.output(print(summary_table)),
    "",
    "Why window-to-gene aggregation was used:",
    "The 25 bp window count strategy first quantifies insertion signals across local genomic windows.",
    "The windows overlapping CDS regions are then aggregated to gene level.",
    "This provides gene-level interpretation while still being closer to the insertion-site logic of Tn-seq than direct HTSeq-style feature counting."
  ),
  con = file.path(OUT_DIR, "window_gene_Tnseq_DESeq2_result_summary.txt")
)

cat("\nDone: window-to-gene Tn-seq DESeq2 analysis finished.\n")
cat("Results saved to:", OUT_DIR, "\n")
