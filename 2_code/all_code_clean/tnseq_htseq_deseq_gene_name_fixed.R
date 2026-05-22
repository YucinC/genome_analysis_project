# ============================================================
# Tn-seq DESeq2 analysis using HTSeq-count output
# Reference: Canu assembly
# Annotation: Prokka TSV from the same Canu assembly
#
# Main fix compared with the previous Tn-seq script:
# 1. Use the same annotation strategy as the RNA-seq script:
#      HTSeq feature/locus_tag -> Prokka locus_tag
# 2. Keep locus_tag only as the internal matching key.
# 3. Use gene_name + product/function for plot labels.
# 4. Add strict annotation-overlap diagnostics so failed annotation matching
#    cannot silently produce unnamed_gene labels.
# ============================================================

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

# HTSeq-count files should be under:
#   COUNT_ROOT/Tn-Seq_BHI
#   COUNT_ROOT/Tn-Seq_HSerum
#   COUNT_ROOT/Tn-Seq_Serum
COUNT_ROOT <- "../tnseq_counting/htseq_canu_bowtie2"

# IMPORTANT:
# This must be the Prokka TSV generated from the SAME Canu annotation/GFF
# that was used for HTSeq-count.
ANNOTATION_FILE <- "../E745_canu.tsv"

OUT_DIR <- "deseq_htseq_trimmed_results_canu_prokka_fixed"
PLOT_DIR <- file.path(OUT_DIR, "plots")
DEBUG_DIR <- file.path(OUT_DIR, "annotation_debug")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DEBUG_DIR, recursive = TRUE, showWarnings = FALSE)

conditions <- c("Tn-Seq_BHI", "Tn-Seq_HSerum", "Tn-Seq_Serum")

LFC_CUTOFF <- 1
PADJ_CUTOFF <- 0.05
N_LABEL_EACH_DIRECTION <- 5

cat("Count root:\n", COUNT_ROOT, "\n")
cat("Annotation file:\n", ANNOTATION_FILE, "\n")
cat("Output directory:\n", OUT_DIR, "\n\n")

if (!file.exists(ANNOTATION_FILE)) {
  stop("Annotation file not found. Please check ANNOTATION_FILE.")
}

# -----------------------------
# 2. Helper functions
# -----------------------------

clean_htseq_gene_id <- function(x) {
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
  x <- sub("\\.canu\\.bowtie2\\.htseq\\.counts\\.txt$", "", x)
  x <- sub("\\.canu\\.htseq\\.counts\\.txt$", "", x)
  x <- sub("\\.htseq\\.counts\\.txt$", "", x)
  x <- sub("\\.counts\\.txt$", "", x)
  x <- sub("\\.txt$", "", x)
  x <- sub("_pass$", "", x)
  return(x)
}

shorten_text <- function(x, n = 55) {
  x <- as.character(x)
  ifelse(!is.na(x) & x != "", substr(x, 1, n), NA)
}

safe_pdf_save <- function(filename, plot, width, height) {
  try(
    ggsave(filename = filename, plot = plot, width = width, height = height),
    silent = TRUE
  )
}

comparison_label <- function(comparison_name) {
  x <- gsub("_vs_", " / ", comparison_name)
  x <- gsub("Tn-Seq_", "", x)
  return(x)
}

# -----------------------------
# 3. Find HTSeq-count files and define metadata
# -----------------------------

count_files <- unlist(lapply(conditions, function(cond) {
  list.files(
    file.path(COUNT_ROOT, cond),
    pattern = "\\.counts\\.txt$",
    full.names = TRUE,
    recursive = TRUE
  )
}))

count_files <- sort(count_files)

if (length(count_files) == 0) {
  stop("No HTSeq count files found. Please check COUNT_ROOT and file pattern.")
}

sample_names <- clean_sample_name(count_files)
dataset <- sub(paste0(".*(", paste(conditions, collapse = "|"), ").*"), "\\1", count_files)

metadata <- data.frame(
  sample = sample_names,
  dataset = dataset,
  condition = dataset,
  file = count_files,
  stringsAsFactors = FALSE
)

metadata$condition <- factor(metadata$condition, levels = conditions)

if (any(is.na(metadata$condition))) {
  stop("Some samples could not be assigned to Tn-seq conditions. Please check folder names.")
}

rownames(metadata) <- metadata$sample

cat("Sample metadata:\n")
print(metadata)

write.csv(
  metadata,
  file.path(OUT_DIR, "sample_metadata.csv"),
  row.names = TRUE,
  quote = FALSE
)

# -----------------------------
# 4. Read HTSeq-count files
# -----------------------------

read_htseq_count_with_summary <- function(file) {
  df <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("gene_id", "count")
  )

  df$gene_id_original <- df$gene_id
  df$gene_id <- clean_htseq_gene_id(df$gene_id)
  df$count <- as.numeric(df$count)
  df$count[is.na(df$count)] <- 0

  assigned <- df[!grepl("^__", df$gene_id), c("gene_id", "count")]
  summary_rows <- df[grepl("^__", df$gene_id), ]

  get_summary_value <- function(name) {
    v <- summary_rows$count[summary_rows$gene_id == name]
    if (length(v) == 0) return(0)
    return(v[1])
  }

  total_in_count_file <- sum(df$count)
  assigned_to_features <- sum(assigned$count)

  summary_df <- data.frame(
    sample = clean_sample_name(file),
    count_file = file,
    total_in_count_file = total_in_count_file,
    assigned_to_features = assigned_to_features,
    no_feature = get_summary_value("__no_feature"),
    ambiguous = get_summary_value("__ambiguous"),
    too_low_aQual = get_summary_value("__too_low_aQual"),
    not_aligned = get_summary_value("__not_aligned"),
    alignment_not_unique = get_summary_value("__alignment_not_unique"),
    assigned_rate_in_count_file_percent = round(100 * assigned_to_features / total_in_count_file, 2),
    stringsAsFactors = FALSE
  )

  return(list(counts = assigned, summary = summary_df))
}

htseq_list <- lapply(metadata$file, read_htseq_count_with_summary)
count_list <- lapply(htseq_list, function(x) x$counts)
names(count_list) <- metadata$sample

htseq_summary <- bind_rows(lapply(htseq_list, function(x) x$summary)) %>%
  left_join(metadata[, c("sample", "condition")], by = "sample") %>%
  select(condition, everything())

write.csv(
  htseq_summary,
  file.path(OUT_DIR, "htseq_counting_summary_from_count_files.csv"),
  row.names = FALSE,
  quote = FALSE
)

all_genes <- Reduce(union, lapply(count_list, function(x) x$gene_id))

cat("\nNumber of HTSeq features detected across all count files:", length(all_genes), "\n")

count_matrix <- matrix(
  0L,
  nrow = length(all_genes),
  ncol = length(count_list),
  dimnames = list(all_genes, names(count_list))
)

for (sample in names(count_list)) {
  df <- count_list[[sample]]
  count_matrix[df$gene_id, sample] <- round(df$count)
}

count_matrix <- count_matrix[, rownames(metadata)]

write.csv(
  as.data.frame(count_matrix),
  file.path(OUT_DIR, "raw_count_matrix_locus_internal.csv"),
  quote = FALSE
)

cat("\nRaw count matrix dimension:\n")
print(dim(count_matrix))
cat("\nLibrary sizes based on assigned Tn-seq insertion counts:\n")
print(colSums(count_matrix))

# -----------------------------
# 5. Read and process Prokka annotation
# -----------------------------

annotation_raw <- read.delim(
  ANNOTATION_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

cat("\nAnnotation table columns:\n")
print(colnames(annotation_raw))

if (!("locus_tag" %in% colnames(annotation_raw))) {
  stop("The annotation file does not contain a 'locus_tag' column. Please check the Prokka TSV format.")
}
if (!("gene" %in% colnames(annotation_raw))) annotation_raw$gene <- NA
if (!("product" %in% colnames(annotation_raw))) annotation_raw$product <- NA
if (!("ftype" %in% colnames(annotation_raw))) annotation_raw$ftype <- NA
if (!("COG" %in% colnames(annotation_raw))) annotation_raw$COG <- NA
if (!("EC_number" %in% colnames(annotation_raw))) annotation_raw$EC_number <- NA

annotation_df <- annotation_raw %>%
  mutate(
    gene_id = clean_htseq_gene_id(locus_tag),
    gene_name = ifelse(!is.na(gene) & gene != "" & gene != "-", gene, NA),
    product_name = ifelse(!is.na(product) & product != "" & product != "-", product, NA),
    product_short = shorten_text(product_name, 55),
    has_gene_name = !is.na(gene_name) & gene_name != "",
    has_product_name = !is.na(product_name) & product_name != "",
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
  arrange(gene_id, annotation_priority) %>%
  group_by(gene_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    plot_label = case_when(
      has_gene_name & has_product_name ~ paste0(gene_name, " / ", product_short),
      has_gene_name ~ gene_name,
      has_product_name ~ product_short,
      TRUE ~ gene_id
    )
  ) %>%
  select(
    gene_id,
    locus_tag,
    ftype,
    gene_name,
    product_name,
    product_short,
    plot_label,
    COG,
    EC_number
  )

write.csv(
  annotation_df,
  file.path(OUT_DIR, "prokka_annotation_simplified_canu.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Annotation matching diagnostics.
count_gene_ids <- rownames(count_matrix)
annotation_gene_ids <- annotation_df$gene_id
overlap_gene_ids <- intersect(count_gene_ids, annotation_gene_ids)
missing_from_annotation <- setdiff(count_gene_ids, annotation_gene_ids)

match_summary <- data.frame(
  count_features = length(count_gene_ids),
  annotation_features = length(annotation_gene_ids),
  matched_features = length(overlap_gene_ids),
  unmatched_count_features = length(missing_from_annotation),
  match_rate_percent = round(100 * length(overlap_gene_ids) / length(count_gene_ids), 2),
  stringsAsFactors = FALSE
)

write.csv(
  match_summary,
  file.path(DEBUG_DIR, "annotation_matching_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(first_50_count_gene_ids = head(count_gene_ids, 50)),
  file.path(DEBUG_DIR, "first_50_htseq_gene_ids.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(first_50_annotation_gene_ids = head(annotation_gene_ids, 50)),
  file.path(DEBUG_DIR, "first_50_annotation_gene_ids.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(unmatched_count_gene_id = head(missing_from_annotation, 200)),
  file.path(DEBUG_DIR, "unmatched_htseq_gene_ids_first200.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nAnnotation matching summary:\n")
print(match_summary)

if (length(overlap_gene_ids) / length(count_gene_ids) < 0.2) {
  stop(
    paste0(
      "Annotation matching rate is very low (",
      match_summary$match_rate_percent,
      "%). This usually means HTSeq-count and Prokka TSV are not from the same GFF/assembly, ",
      "or the feature ID format is different. Check files in: ", DEBUG_DIR
    )
  )
}

# A display-name version of count matrix for human reading.
count_annotation_for_matrix <- data.frame(gene_id = rownames(count_matrix), stringsAsFactors = FALSE) %>%
  left_join(annotation_df, by = "gene_id") %>%
  mutate(
    display_id = case_when(
      !is.na(gene_name) & gene_name != "" ~ gene_name,
      !is.na(product_short) & product_short != "" ~ product_short,
      TRUE ~ gene_id
    )
  )

count_matrix_display <- as.data.frame(count_matrix)
count_matrix_display <- cbind(
  gene_id = count_annotation_for_matrix$display_id,
  product_name = ifelse(is.na(count_annotation_for_matrix$product_name), "", count_annotation_for_matrix$product_name),
  locus_tag_internal = rownames(count_matrix),
  count_matrix_display
)

write.csv(
  count_matrix_display,
  file.path(OUT_DIR, "raw_count_matrix_gene_name_annotated.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# 6. Prefilter and run DESeq2
# -----------------------------

keep <- rowSums(count_matrix) >= 10
count_matrix_filtered <- count_matrix[keep, ]

cat("\nFeatures before filtering:", nrow(count_matrix), "\n")
cat("Features after filtering:", nrow(count_matrix_filtered), "\n")

metadata$condition <- factor(metadata$condition, levels = conditions)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix_filtered),
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

saveRDS(dds, file.path(OUT_DIR, "dds_htseq_trimmed_tnseq_canu.rds"))

size_factor_table <- data.frame(
  sample = colnames(count_matrix_filtered),
  condition = metadata$condition,
  raw_assigned_counts_after_filtering = colSums(count_matrix_filtered),
  size_factor = sizeFactors(dds),
  stringsAsFactors = FALSE
)

write.csv(
  size_factor_table,
  file.path(OUT_DIR, "deseq2_size_factors.csv"),
  row.names = FALSE,
  quote = FALSE
)

normalized_counts <- counts(dds, normalized = TRUE)

normalized_annotation_for_matrix <- data.frame(gene_id = rownames(normalized_counts), stringsAsFactors = FALSE) %>%
  left_join(annotation_df, by = "gene_id") %>%
  mutate(
    display_id = case_when(
      !is.na(gene_name) & gene_name != "" ~ gene_name,
      !is.na(product_short) & product_short != "" ~ product_short,
      TRUE ~ gene_id
    )
  )

normalized_counts_display <- as.data.frame(normalized_counts)
normalized_counts_display <- cbind(
  gene_id = normalized_annotation_for_matrix$display_id,
  product_name = ifelse(is.na(normalized_annotation_for_matrix$product_name), "", normalized_annotation_for_matrix$product_name),
  locus_tag_internal = rownames(normalized_counts),
  normalized_counts_display
)

write.csv(
  normalized_counts_display,
  file.path(OUT_DIR, "normalized_counts_gene_name_annotated.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# 7. Global QC plots
# -----------------------------

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(OUT_DIR, "vst_htseq_trimmed_tnseq_canu.rds"))

pca_matrix <- t(assay(vsd))
pca <- prcomp(pca_matrix)
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = metadata[rownames(pca$x), "condition"],
  stringsAsFactors = FALSE
)

write.csv(pca_df, file.path(OUT_DIR, "PCA_coordinates.csv"), row.names = FALSE, quote = FALSE)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3.5, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3) +
  labs(
    title = "PCA plot of Tn-seq HTSeq-count samples",
    subtitle = "VST-transformed insertion counts",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right") +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))

ggsave(file.path(PLOT_DIR, "PCA_plot_HTSeq_trimmed_Tnseq_labeled.png"), p_pca, width = 7.5, height = 7, dpi = 300)
safe_pdf_save(file.path(PLOT_DIR, "PCA_plot_HTSeq_trimmed_Tnseq_labeled.pdf"), p_pca, width = 7.5, height = 7)

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

annotation_col <- data.frame(condition = colData(vsd)$condition)
rownames(annotation_col) <- colnames(vsd)

write.csv(sample_dist_matrix, file.path(OUT_DIR, "sample_distance_matrix.csv"), quote = FALSE)

png(file.path(PLOT_DIR, "Sample_distance_heatmap_HTSeq_trimmed_Tnseq.png"), width = 2100, height = 1800, res = 300)
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
  pdf(file.path(PLOT_DIR, "Sample_distance_heatmap_HTSeq_trimmed_Tnseq.pdf"), width = 7, height = 6)
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
# 8. Pairwise comparisons and plots
# -----------------------------

comparisons <- list(
  Serum_vs_BHI = c("condition", "Tn-Seq_Serum", "Tn-Seq_BHI"),
  HSerum_vs_BHI = c("condition", "Tn-Seq_HSerum", "Tn-Seq_BHI"),
  Serum_vs_HSerum = c("condition", "Tn-Seq_Serum", "Tn-Seq_HSerum")
)

annotate_result <- function(res_df, annotation_df) {
  res_df$gene_id <- clean_htseq_gene_id(res_df$gene_id)

  res_annotated <- res_df %>%
    left_join(annotation_df, by = "gene_id") %>%
    mutate(
      gene_name = ifelse(is.na(gene_name), "", gene_name),
      product_name = ifelse(is.na(product_name), "", product_name),
      product_short = ifelse(is.na(product_short), "", product_short),
      plot_label = case_when(
        gene_name != "" & product_short != "" ~ paste0(gene_name, " / ", product_short),
        gene_name != "" ~ gene_name,
        product_short != "" ~ product_short,
        TRUE ~ gene_id
      ),
      display_gene_id = case_when(
        gene_name != "" ~ gene_name,
        product_short != "" ~ product_short,
        TRUE ~ gene_id
      ),
      annotation_status = case_when(
        gene_name != "" & product_name != "" ~ "gene_and_product",
        gene_name != "" ~ "gene_only",
        product_name != "" ~ "product_only",
        TRUE ~ "not_matched_or_no_annotation"
      )
    ) %>%
    relocate(display_gene_id, gene_name, product_name, product_short, plot_label, .after = gene_id)

  n_total <- nrow(res_annotated)
  n_matched <- sum(res_annotated$annotation_status != "not_matched_or_no_annotation", na.rm = TRUE)
  message("Annotation matching summary for result table:")
  message("  Total result features: ", n_total)
  message("  Features with gene/product annotation: ", n_matched)
  message("  Annotation rate: ", round(100 * n_matched / n_total, 2), "%")

  return(res_annotated)
}

classify_tnseq_change <- function(df, lfc_cutoff = 1, padj_cutoff = 0.05) {
  df$tnseq_change <- "Not significant"

  higher_idx <- !is.na(df$padj) & df$padj < padj_cutoff &
    !is.na(df$log2FoldChange) & df$log2FoldChange >= lfc_cutoff

  lower_idx <- !is.na(df$padj) & df$padj < padj_cutoff &
    !is.na(df$log2FoldChange) & df$log2FoldChange <= -lfc_cutoff

  df$tnseq_change[higher_idx] <- "Higher mutant abundance"
  df$tnseq_change[lower_idx] <- "Lower mutant abundance"

  df$tnseq_change <- factor(
    df$tnseq_change,
    levels = c("Lower mutant abundance", "Not significant", "Higher mutant abundance")
  )

  return(df)
}

select_label_genes <- function(df, n_each_direction = 5, lfc_cutoff = 1, padj_cutoff = 0.05) {
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
    distinct(gene_id, .keep_all = TRUE)
}

summary_rows <- list()

for (name in names(comparisons)) {
  cat("\nRunning comparison:", name, "\n")

  res <- results(dds, contrast = comparisons[[name]])
  res <- res[order(res$padj), ]

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

  res_annotated <- annotate_result(res_df, annotation_df)
  res_annotated <- classify_tnseq_change(res_annotated, LFC_CUTOFF, PADJ_CUTOFF)

  # Human-facing table: gene display columns first, internal locus kept later for traceability.
  res_export <- res_annotated %>%
    transmute(
      gene_id = display_gene_id,
      gene_name = gene_name,
      product_name = product_name,
      product_short = product_short,
      plot_label = plot_label,
      baseMean = baseMean,
      log2FoldChange = log2FoldChange,
      lfcSE = lfcSE,
      stat = stat,
      pvalue = pvalue,
      padj = padj,
      tnseq_change = tnseq_change,
      annotation_status = annotation_status,
      locus_tag_internal = gene_id,
      ftype = ftype,
      COG = COG,
      EC_number = EC_number
    )

  annotation_check <- res_export %>% count(annotation_status)
  write.csv(annotation_check, file.path(DEBUG_DIR, paste0(name, "_annotation_status_summary.csv")), row.names = FALSE, quote = FALSE)

  write.csv(
    res_export,
    file.path(OUT_DIR, paste0(name, "_HTSeq_trimmed_DESeq2_results_annotated.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  sig_padj <- res_export %>% filter(!is.na(padj), padj < PADJ_CUTOFF)
  sig_lfc <- res_export %>% filter(!is.na(padj), padj < PADJ_CUTOFF, !is.na(log2FoldChange), abs(log2FoldChange) >= LFC_CUTOFF)

  write.csv(sig_padj, file.path(OUT_DIR, paste0(name, "_HTSeq_trimmed_significant_padj0.05_annotated.csv")), row.names = FALSE, quote = FALSE)
  write.csv(sig_lfc, file.path(OUT_DIR, paste0(name, "_HTSeq_trimmed_significant_padj0.05_log2FC1_annotated.csv")), row.names = FALSE, quote = FALSE)

  label_genes <- select_label_genes(res_annotated, N_LABEL_EACH_DIRECTION, LFC_CUTOFF, PADJ_CUTOFF)
  label_genes_export <- label_genes %>%
    transmute(
      gene_id = display_gene_id,
      gene_name = gene_name,
      product_name = product_name,
      product_short = product_short,
      plot_label = plot_label,
      baseMean = baseMean,
      log2FoldChange = log2FoldChange,
      padj = padj,
      tnseq_change = tnseq_change,
      locus_tag_internal = gene_id,
      annotation_status = annotation_status
    )

  write.csv(label_genes_export, file.path(OUT_DIR, paste0(name, "_HTSeq_trimmed_labeled_genes_in_plots.csv")), row.names = FALSE, quote = FALSE)

  n_higher <- sum(res_annotated$tnseq_change == "Higher mutant abundance", na.rm = TRUE)
  n_lower <- sum(res_annotated$tnseq_change == "Lower mutant abundance", na.rm = TRUE)

  ma_df <- res_annotated %>%
    mutate(
      baseMean_for_plot = ifelse(baseMean <= 0 | is.na(baseMean), NA, baseMean),
      label_this_gene = gene_id %in% label_genes$gene_id
    ) %>%
    filter(!is.na(baseMean_for_plot), !is.na(log2FoldChange))

  p_ma <- ggplot(ma_df, aes(x = baseMean_for_plot, y = log2FoldChange, color = tnseq_change)) +
    geom_point(alpha = 0.7, size = 1.4) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_label_repel(
      data = ma_df %>% filter(label_this_gene),
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
      subtitle = paste0("Labels use Prokka gene name + product/function. Threshold: padj < ", PADJ_CUTOFF, " and |log2FC| >= ", LFC_CUTOFF, "\n",
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

  ggsave(file.path(PLOT_DIR, paste0(name, "_HTSeq_trimmed_MA_plot_annotated.png")), p_ma, width = 8.5, height = 6, dpi = 300)
  safe_pdf_save(file.path(PLOT_DIR, paste0(name, "_HTSeq_trimmed_MA_plot_annotated.pdf")), p_ma, width = 8.5, height = 6)

  volcano_df <- res_annotated %>%
    mutate(
      padj_for_plot = case_when(is.na(padj) ~ NA_real_, padj == 0 ~ .Machine$double.xmin, TRUE ~ padj),
      neg_log10_padj = -log10(padj_for_plot),
      label_this_gene = gene_id %in% label_genes$gene_id
    ) %>%
    filter(!is.na(log2FoldChange), !is.na(neg_log10_padj))

  p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = tnseq_change)) +
    geom_point(alpha = 0.75, size = 1.4) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_label_repel(
      data = volcano_df %>% filter(label_this_gene),
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
                        "Labels use Prokka gene name + product/function. Threshold: padj < ", PADJ_CUTOFF, " and |log2FC| >= ", LFC_CUTOFF),
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

  ggsave(file.path(PLOT_DIR, paste0(name, "_HTSeq_trimmed_volcano_plot_annotated.png")), p_volcano, width = 8.5, height = 6, dpi = 300)
  safe_pdf_save(file.path(PLOT_DIR, paste0(name, "_HTSeq_trimmed_volcano_plot_annotated.pdf")), p_volcano, width = 8.5, height = 6)

  top20 <- res_export %>% filter(!is.na(padj)) %>% arrange(padj) %>% slice_head(n = 20)
  write.csv(top20, file.path(OUT_DIR, paste0(name, "_HTSeq_trimmed_top20_by_padj_annotated.csv")), row.names = FALSE, quote = FALSE)

  summary_rows[[name]] <- data.frame(
    analysis = "HTSeq_trimmed_gene_level",
    comparison = name,
    tested_features = nrow(res_export),
    annotated_features = sum(res_export$annotation_status != "not_matched_or_no_annotation", na.rm = TRUE),
    significant_padj_0.05 = nrow(sig_padj),
    significant_padj_0.05_log2FC_1 = nrow(sig_lfc),
    higher_mutant_abundance = n_higher,
    lower_mutant_abundance = n_lower,
    labeled_genes = nrow(label_genes_export),
    top_gene_id = ifelse(nrow(top20) > 0, top20$gene_id[1], NA),
    top_plot_label = ifelse(nrow(top20) > 0, top20$plot_label[1], NA),
    top_locus_tag_internal = ifelse(nrow(top20) > 0, top20$locus_tag_internal[1], NA),
    top_log2FC = ifelse(nrow(top20) > 0, top20$log2FoldChange[1], NA),
    top_padj = ifelse(nrow(top20) > 0, top20$padj[1], NA),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 9. Final summary
# -----------------------------

summary_table <- bind_rows(summary_rows)

write.csv(
  summary_table,
  file.path(OUT_DIR, "HTSeq_trimmed_Tnseq_DESeq2_result_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  c(
    "Tn-seq DESeq2 result summary: trimmed HTSeq-count gene-level analysis",
    "=====================================================================",
    "",
    "Interpretation:",
    "This is Tn-seq insertion-count analysis, not RNA-seq expression analysis.",
    "Positive log2FoldChange means higher mutant/insertion abundance in the numerator condition.",
    "Negative log2FoldChange means lower mutant/insertion abundance in the numerator condition.",
    "A significantly lower mutant abundance may suggest that disruption of the corresponding gene reduces fitness in the numerator condition.",
    "",
    paste0("Features before filtering: ", nrow(count_matrix)),
    paste0("Features after filtering: ", nrow(count_matrix_filtered)),
    paste0("Annotation match rate before filtering: ", match_summary$match_rate_percent, "%"),
    "",
    "Important annotation note:",
    "The internal DESeq2 rownames are still locus_tag IDs because they are stable and unique.",
    "Human-facing result tables and plot labels use Prokka gene_name/product annotations.",
    "If some labels still show locus_tag, those features were not matched or lack Prokka gene/product annotation.",
    "Check annotation_debug/ for matching diagnostics.",
    "",
    "Pairwise comparison summary:",
    capture.output(print(summary_table))
  ),
  con = file.path(OUT_DIR, "HTSeq_trimmed_Tnseq_DESeq2_result_summary.txt")
)

cat("\nAnalysis finished successfully.\n")
cat("Results saved to:\n", OUT_DIR, "\n")
cat("Main output files use gene_name/product annotations.\n")
cat("Check annotation diagnostics in:\n", DEBUG_DIR, "\n")
