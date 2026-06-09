#!/usr/bin/env Rscript
# =============================================================================
# cn_plot_cluster.R
# Copy-number genome plot with integrated hierarchical clustering
#
# Usage:
#   Rscript cn_plot_cluster.R \
#     --data          /path/to/dataCHORDOMA_dataMASTER.RData \
#     --gencode       /path/to/gencode19_gns_lite.RData \
#     --samples       /path/to/samples_selection.tsv \
#     --metadata      /path/to/metadata_annotation.tsv \
#     --centromeres   /path/to/hg19_centromeres.txt \
#     --outdir        /path/to/output \
#     --seq-type      WGS \
#     --k             3 \
#     --bin-size      1000000 \
#     --na-thresh     0.2 \
#     --width         20 \
#     --height        15
#Input
# sample_id   ← SampleID / sample / Tumor_Sample_Barcode
# chr         ← chrom / chromosome / seqnames
# start/end   ← startpos / chromStart / endpos
# TCN         ← total_cn / Copy_Number / nMajor+nMinor (derived)
## Mode A — RData (same as before, --samples now optional)
# Rscript cn_plot_cluster.R \
#   --data        dataCHORDOMA.RData \
#   --metadata    oncoprint_meta_annotation.tsv \
#   --centromeres hg19_centromeres.txt \
#   --outdir      results/ --seq-type WGS --k 3

# # Mode B — ASCAT flat output, no sample sheet needed
# Rscript cn_plot_cluster.R \
#   --data        ascat_segments_all_samples.tsv \
#   --metadata    sample_metadata.tsv \
#   --centromeres hg19_centromeres.txt \
#   --outdir      results/ --k 4

# # Mode B — DNA only, subset to specific samples
# Rscript cn_plot_cluster.R \
#   --data        purple_cn_segments.tsv \
#   --metadata    sample_metadata.tsv \
#   --samples     samples_of_interest.tsv \
#   --centromeres hg19_centromeres.txt \
#   --outdir      results/
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(cowplot)
  library(RColorBrewer)
  library(ggplot2)
  library(MultiAssayExperiment)
  library(GenomicRanges)
  library(pheatmap)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("--data"),        type = "character", help = "Path to dataCHORDOMA_dataMASTER.RData [required]"),
  make_option(c("--gencode"),     type = "character", help = "Path to gencode19_gns_lite.RData [required]"),
  make_option(c("--samples"),     type = "character", help = "Path to samples_selection.tsv [required]"),
  make_option(c("--metadata"),    type = "character", help = "Path to metadata_annotation.tsv [required]"),
  make_option(c("--centromeres"), type = "character", help = "Path to hg19_centromeres.txt [required]"),
  make_option(c("--outdir"),      type = "character", default = ".",
              help = "Output directory [default: current dir]"),
  make_option(c("--seq-type"),    type = "character", default = "WGS",
              help = "Sequencing type to plot: WGS, WES, or ALL [default: WGS]"),
  make_option(c("-k", "--k"),     type = "integer",   default = 3L,
              help = "Number of clusters for cutree [default: 3]"),
  make_option(c("--bin-size"),    type = "double",    default = 1e6,
              help = "Bin size in bp for CN matrix [default: 1000000]"),
  make_option(c("--na-thresh"),   type = "double",    default = 0.2,
              help = "Max fraction of NAs per bin before removal [default: 0.2]"),
  make_option(c("--width"),       type = "double",    default = 20,
              help = "Plot width in inches [default: 20]"),
  make_option(c("--height"),      type = "double",    default = 15,
              help = "Plot height in inches [default: 15]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── Validate required args ────────────────────────────────────────────────────
for (arg in c("data", "gencode", "samples", "metadata", "centromeres")) {
  if (is.null(opt[[arg]])) stop("Missing required argument: --", arg)
  if (!file.exists(opt[[arg]])) stop("File not found: ", opt[[arg]])
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

seq_filter <- if (toupper(opt$`seq-type`) == "ALL") NULL else opt$`seq-type`

message("\n=== CN Plot + Clustering ===")
message("  seq-type  : ", ifelse(is.null(seq_filter), "ALL", seq_filter))
message("  k         : ", opt$k)
message("  bin-size  : ", format(opt$`bin-size`, big.mark = ","), " bp")
message("  na-thresh : ", opt$`na-thresh`)
message("  outdir    : ", opt$outdir)

# =============================================================================
# STEP 1 — Load data
# =============================================================================
message("\n[1/5] Loading data...")

load(opt$gencode)
load(opt$data)

colData_df <- as.data.frame(colData(dataCHORDOMA))

sample_selection <- read_tsv(opt$samples, show_col_types = FALSE) %>%
  filter(SelectedSample)

cnvSmpls <- sample_selection %>% filter(TakeCNV)
smplsCN  <- cnvSmpls$ID

cohortCNA <- dataCHORDOMA[, smplsCN, "cnv"]
cohortCNA <- as(cohortCNA[["cnv"]], "GRangesList")

cohort <- as.data.frame(cohortCNA)
cohort$ploidy         <- sample_selection$Ploidy[match(cohort$group_name, sample_selection$ID)]
cohort$purity         <- sample_selection$TumorCellContent[match(cohort$group_name, sample_selection$ID)]
cohort$predicted_sex  <- sample_selection$Sex[match(cohort$group_name, sample_selection$ID)]
cohort$PID            <- sample_selection$PatientID[match(cohort$group_name, sample_selection$ID)]

metadata <- read_tsv(opt$metadata, show_col_types = FALSE)
cohort$localization   <- metadata$LocalizationBasket[match(cohort$PID, metadata$PatientID)]

cohort$ploidy_round <- round(cohort$ploidy)

cohort <- cohort %>%
  mutate(seq_type = case_when(
    grepl("^WES", group_name) ~ "WES",
    grepl("^WGS", group_name) ~ "WGS",
    TRUE ~ "other"
  ))

if (!is.null(seq_filter)) {
  cohort <- cohort %>% filter(seq_type == seq_filter)
  if (nrow(cohort) == 0) stop("No samples remain after seq-type filter: ", seq_filter)
}

message("  Samples loaded: ", length(unique(cohort$group_name)))

# =============================================================================
# STEP 2 — Prepare segment data
# =============================================================================
message("\n[2/5] Preparing segment data...")

x <- cohort %>%
  mutate(cn_diff = TCN - ploidy_round) %>%
  mutate(seqnames = ifelse(seqnames == "X", "23", seqnames)) %>%
  filter(!seqnames %in% c("Y", "24"))

x.m <- x %>%
  mutate(segment_id = seq_len(n())) %>%
  pivot_longer(cols = c(start, end),
               names_to  = "start_end",
               values_to = "position") %>%
  mutate(
    position_mb          = position / 1e6,
    chr                  = factor(paste0("chr", seqnames), levels = paste0("chr", 1:23)),
    log2_fold_change     = log2(TCN / ploidy_round),
    log2_fold_change     = ifelse(is.nan(log2_fold_change), -5, log2_fold_change),
    log2_fold_change_cut = case_when(
      log2_fold_change >  2 ~  2,
      log2_fold_change < -2 ~ -2,
      TRUE                  ~ log2_fold_change
    )
  )

# annotation data frame (one row per sample)
y <- cohort %>%
  dplyr::select(group_name, predicted_sex, purity, ploidy, seq_type, localization) %>%
  distinct() %>%
  mutate(chr = factor("chr1", levels = paste0("chr", 1:23)))

# =============================================================================
# STEP 3 — Centromeres
# =============================================================================
message("\n[3/5] Loading centromeres from: ", opt$centromeres)

centromeres <- read.table(opt$centromeres, header = TRUE, comment.char = "") %>%
  dplyr::rename(chr = chrom) %>%
  filter(chr != "chrY") %>%
  mutate(
    chr = gsub("chrX", "chr23", chr),
    chr = factor(chr, levels = paste0("chr", 1:23))
  ) %>%
  dplyr::select(chr, chromStart, chromEnd) %>%
  pivot_longer(cols = c(chromStart, chromEnd),
               names_to  = "start_end",
               values_to = "position_mb") %>%
  mutate(position_mb = position_mb / 1e6)

# colour scale
mapal <- colorRampPalette(brewer.pal(11, "RdBu"))(256)

# =============================================================================
# STEP 4 — Build CN matrix + cluster (used by the plot function)
# =============================================================================
build_cn_matrix <- function(x_sub, bin_size, na_thresh) {

  genome_bins <- x_sub %>%
    group_by(seqnames) %>%
    summarise(max_pos = max(end), .groups = "drop") %>%
    rowwise() %>%
    mutate(bin_start = list(seq(1, max_pos, by = bin_size))) %>%
    unnest(bin_start) %>%
    mutate(
      bin_end = bin_start + bin_size - 1,
      bin_id  = paste0(seqnames, ":", bin_start)
    )

  cn_matrix <- x_sub %>%
    mutate(log2fc = log2(pmax(TCN, 0.01) / ploidy_round)) %>%
    group_by(group_name) %>%
    group_modify(~ {
      sample_seg <- .x
      genome_bins %>%
        rowwise() %>%
        mutate(log2fc = {
          overlap <- sample_seg %>%
            filter(seqnames == seqnames,
                   start    <= bin_end,
                   end      >= bin_start)
          if (nrow(overlap) == 0) NA_real_ else mean(overlap$log2fc)
        }) %>%
        ungroup()
    }) %>%
    ungroup()

  cn_wide <- cn_matrix %>%
    dplyr::select(group_name, bin_id, log2fc) %>%
    pivot_wider(names_from = bin_id, values_from = log2fc) %>%
    column_to_rownames("group_name")

  cn_wide <- cn_wide[, colMeans(is.na(cn_wide)) < na_thresh]
  cn_wide <- apply(cn_wide, 2, function(col) {
    col[is.na(col)] <- median(col, na.rm = TRUE)
    col
  }) %>% as.data.frame()

  cn_wide
}

# =============================================================================
# STEP 5 — Plot function (clustering integrated)
# =============================================================================
make_full_plot_ordered <- function(seq        = NULL,
                                   bin_size   = 1e6,
                                   k          = 3L,
                                   na_thresh  = 0.2) {

  # ── Filter ─────────────────────────────────────────────────────────────
  y_sub   <- y   %>% { if (!is.null(seq)) filter(., seq_type == seq) else . }
  x_sub   <- x   %>% { if (!is.null(seq)) filter(., seq_type == seq) else . }
  x.m_sub <- x.m %>% { if (!is.null(seq)) filter(., seq_type == seq) else . }

  message("  Building CN matrix for clustering (", nrow(x_sub), " segments)...")

  # ── Cluster ─────────────────────────────────────────────────────────────
  cn_wide   <- build_cn_matrix(x_sub, bin_size, na_thresh)
  dist_mat  <- dist(cn_wide, method = "euclidean")
  hclust_cn <- hclust(dist_mat, method = "ward.D2")

  sample_order <- hclust_cn$labels[hclust_cn$order]
  clusters     <- cutree(hclust_cn, k = k)
  cluster_df   <- data.frame(
    group_name = names(clusters),
    cluster    = factor(clusters)
  )

  cluster_colors <- setNames(
    RColorBrewer::brewer.pal(max(3L, k), "Set2")[seq_len(k)],
    levels(cluster_df$cluster)
  )

  message("  Cluster sizes: ",
          paste(table(cluster_df$cluster), collapse = " | "))

  # ── Apply order ──────────────────────────────────────────────────────────
  y_sub <- y_sub %>%
    left_join(cluster_df, by = "group_name") %>%
    mutate(group_name = factor(group_name, levels = sample_order))

  x.m_sub <- x.m_sub %>%
    left_join(cluster_df, by = "group_name") %>%
    mutate(group_name = factor(group_name, levels = sample_order))

  blank_data_sub <- x.m_sub %>%
    group_by(chr) %>%
    summarise(position_mb = ceiling(max(position_mb)), .groups = "drop")

  # ── Blank facet template ─────────────────────────────────────────────────
  p_blank <- ggplot(y_sub) +
    facet_grid(. ~ group_name, space = "free", scales = "free") +
    theme_cowplot() +
    panel_border(colour = "gray80", size = 0.5, linetype = 1) +
    theme(
      strip.placement  = "outside",
      strip.text.x     = element_text(angle = 90, size = 1),
      axis.title       = element_blank(),
      axis.text.y      = element_text(size = 7),
      axis.text.x      = element_blank(),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      panel.spacing    = unit(0.1, "lines"),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))

  # ── Annotation tiles ─────────────────────────────────────────────────────
  p_cluster <- p_blank +
    geom_tile(aes("1", "Cluster", fill = cluster)) +
    scale_fill_manual(name = "Cluster", values = cluster_colors)

  p_sex <- p_blank +
    geom_tile(aes("1", "Sex", fill = predicted_sex)) +
    scale_fill_manual(name   = "Sex",
                      values = c("female" = "darkgoldenrod1", "male" = "darkcyan"))

  p_purity <- p_blank +
    geom_tile(aes("1", "Purity", fill = purity)) +
    scale_fill_gradientn(name    = "Purity",
                         colours = brewer.pal(9, "Blues"),
                         limits  = c(0, 1))

  p_ploidy <- p_blank +
    geom_tile(aes("1", "Ploidy", fill = ploidy)) +
    scale_fill_gradientn(name    = "Ploidy",
                         colours = brewer.pal(9, "Greens"))

  p_loc <- p_blank +
    geom_tile(aes("1", "Loc", fill = localization))

  # ── Main CN genome plot ──────────────────────────────────────────────────
  p_cn <- ggplot(x.m_sub) +
    geom_path(aes(x = 0, y = position_mb,
                  group = segment_id,
                  color = log2_fold_change_cut),
              linewidth = 10) +
    geom_blank(data = blank_data_sub,
               aes(x = 0, y = position_mb)) +
    geom_hline(data = centromeres,
               aes(yintercept = position_mb),
               color = "gray80", linewidth = 0.5) +
    facet_grid(chr ~ group_name, space = "free_y", scales = "free_y", switch = "y") +
    scale_color_gradientn(colours = rev(mapal), name = "log2 fold change") +
    scale_y_reverse(
      breaks = function(lim) c(ceiling(min(lim)), floor(max(lim))),
      expand = expansion(mult = 0)
    ) +
    theme_cowplot() +
    panel_border(colour = "gray80", size = 0.5, linetype = 1) +
    theme(
      strip.placement   = "outside",
      strip.text.x      = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      axis.title.x      = element_blank(),
      axis.text.x       = element_blank(),
      axis.line         = element_blank(),
      axis.ticks.x      = element_blank(),
      axis.text.y       = element_text(size = 7),
      panel.spacing     = unit(0.1, "lines")
    ) +
    ylab("Genomic position [Mb]")

  # ── Combine ──────────────────────────────────────────────────────────────
  plot_grid(
    p_cluster, p_sex, p_purity, p_ploidy, p_loc, p_cn,
    rel_heights = c(2, 2, 2, 2, 2, 30),
    ncol        = 1,
    align       = "v",
    axis        = "lr"
  )
}

# =============================================================================
# Generate outputs
# =============================================================================
message("\n[4/5] Generating CN genome plot...")

p_out <- make_full_plot_ordered(
  seq       = seq_filter,
  bin_size  = opt$`bin-size`,
  k         = opt$k,
  na_thresh = opt$`na-thresh`
)

seq_tag   <- ifelse(is.null(seq_filter), "ALL", seq_filter)
cn_pdf    <- file.path(opt$outdir, paste0("cn_log2fc_", seq_tag, ".pdf"))
ggsave(cn_pdf, p_out, width = opt$width, height = opt$height)
message("  Saved: ", cn_pdf)

# ── Heatmap ───────────────────────────────────────────────────────────────────
message("\n[5/5] Generating clustering heatmap...")

x_sub_hm <- x %>% { if (!is.null(seq_filter)) filter(., seq_type == seq_filter) else . }
cn_wide_hm <- build_cn_matrix(x_sub_hm, opt$`bin-size`, opt$`na-thresh`)

dist_mat_hm  <- dist(cn_wide_hm, method = "euclidean")
hclust_hm    <- hclust(dist_mat_hm, method = "ward.D2")
clusters_hm  <- cutree(hclust_hm, k = opt$k)
cluster_df_hm <- data.frame(
  group_name = names(clusters_hm),
  cluster    = factor(clusters_hm)
)

anno_row <- y %>%
  dplyr::select(group_name, seq_type, predicted_sex, purity, ploidy, localization) %>%
  distinct() %>%
  left_join(cluster_df_hm, by = "group_name") %>%
  column_to_rownames("group_name") %>%
  .[rownames(cn_wide_hm), , drop = FALSE]

heatmap_pdf <- file.path(opt$outdir, paste0("cn_clustering_heatmap_", seq_tag, ".pdf"))
pdf(heatmap_pdf, width = 12, height = 8)
pheatmap(
  as.matrix(cn_wide_hm),
  cluster_rows   = hclust_hm,
  cluster_cols   = TRUE,
  annotation_row = anno_row,
  show_colnames  = FALSE,
  show_rownames  = TRUE,
  color          = rev(mapal),
  breaks         = seq(-2, 2, length.out = 257),
  fontsize_row   = 4,
  border_color   = NA
)
dev.off()
message("  Saved: ", heatmap_pdf)

# ── Write cluster assignments ─────────────────────────────────────────────────
cluster_tsv <- file.path(opt$outdir, paste0("cn_cluster_assignments_", seq_tag, ".tsv"))
cluster_df_hm %>%
  rename(sample_id = group_name) %>%
  write_tsv(cluster_tsv)
message("  Saved: ", cluster_tsv)

message("\n[DONE] All outputs written to: ", opt$outdir)