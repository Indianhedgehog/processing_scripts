## the script contains the data for chromothripsis analysis using the OTP data
## The SV data was copied in a  local folder because of the convineince and the 
# to avoid writing unnecessary codes. The old code bloackhave been commented out 
# so that it can used later if required . The CN data is being read from the Robject.
# The ID's are then finally being used for matching the CN and SV names and do the core analysis.

## At present the inversion cannot be subclassified into head and tail because of the way the data is being processed. 
#The strand information is being extracted from the source column in the SV file. This is not ideal but it is a workaround for now. 
#The code can be modified later to extract the strand information from the VCF file if required.

suppressMessages(library(ShatterSeek))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(qs))
suppressMessages(library(readr))
suppressMessages(library(MultiAssayExperiment))
suppressMessages(library(RaggedExperiment))
suppressMessages(library(data.table))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(qs))
suppressMessages(library(ggplot2))

#---------  SV data proicessing ---------#

setwd("/path/to/data/sv_files")

sample_file <- read_tsv("/path/to/data/cna_sv_path.tsv")
sample_file$full_pid <- paste0(sample_file$pid,"_", sample_file$tumor_control)

filenames_sv <- tibble(filenames = list.files(".", pattern = "svs_.*_filtered_dedup_somatic\\.tsv$"))
filenames_sv$pid <- sapply(strsplit(filenames_sv$filenames, "_"), "[", 2)

hipo21 <- filenames_sv %>% 
  filter(!grepl("OE0246", pid))
hipo21$pid <- basename(hipo21$filenames)
hipo21$pid <- gsub("svs_|_filtered_dedup_somatic.tsv", "", hipo21$pid)
hipo21$full_pid <- hipo21$pid
hipo21$full_pid  <- gsub("-", "_", hipo21$full_pid)


uk <- filenames_sv %>% 
  filter(grepl("OE0246", pid))

uk$pid <- basename(uk$filenames)
uk$pid <- gsub("svs_|_filtered_dedup_somatic.tsv", "", uk$pid)
uk$full_pid <- uk$pid
uk$full_pid  <- gsub("-", "_", uk$full_pid)


filenames_sv <- rbind(hipo21 ,uk)
filenames_sv$filenames <- paste0("/path/to/data/sv_files/",filenames_sv$filenames)

sample_selection <- qread("../sample_selection_corrected.qs")
sample_selection$full_pid <- paste0(sample_selection$ProjectID,"_", sample_selection$PatientID,"_",sample_selection$TumorID,"_", sample_selection$ControlID)
sample_selection$pid <- paste0(sample_selection$ProjectID,"_", sample_selection$PatientID)
final_samples<- sample_selection %>% filter(full_pid %in% filenames_sv$full_pid)


#---------- CNA data ---------------#
load("/path/to/data/annotations/gencode19_gns_lite.RData")
load("/path/to/data/objects/dataCHORDOMA_dataMASTER_250828.RData")
colData <- as.data.frame(colData(dataCHORDOMA))
smpls <- rownames(colData)

smplsCN <- final_samples$ID # Only include samples with sufficient CNA data
cohortCNA <- dataCHORDOMA[,smplsCN,"cnv"]
cohortCNA <- as(cohortCNA[["cnv"]], "GRangesList")

cna <- as.data.frame(cohortCNA)
cna <- cna %>% filter(!grepl("WES", group_name))

##filter the sv data based on the copy number samples available
latest_data <- filenames_sv %>% filter(full_pid %in% sample_selection$full_pid)

cna_list <- split(cna, cna$group_name)

latest_data$id <- sample_selection$ID[sample_selection$full_pid %in% latest_data$full_pid]
latest_data <- latest_data[order(latest_data$id), ]
##########################################

#---    Core fucntion of the shatterseek analysis. ------#

valid_chroms <- c(as.character(1:22), "X")
output_dir <- "/path/to/output/shatterseek_results"


for (i in seq_along(latest_data$filenames)) {

  seg.df     <- latest_data$filenames[i]
  basename_i <- tools::file_path_sans_ext(basename(seg.df))
  sv_id      <- latest_data$id[i] # extract shared ID

  print(paste("ID:", sv_id))

  # --- Match CNA by ID, not by position ---
  cna_name <- names(cna_list)[str_detect(names(cna_list), sv_id)]

  if (length(cna_name) == 0) {
    warning(paste("No matching CNA found for:", sv_id, "— skipping"))
    next
  }
  if (length(cna_name) > 1) {
    warning(paste("Multiple CNA matches for:", sv_id, "— taking first"))
    cna_name <- cna_name[1]
  }

  df3 <- cna_list[[cna_name]]
  print(paste("Matched CNA sample:", cna_name))

  # --- Load & filter SV data ---
  df2 <- read_tsv(seg.df)
  df2 <- df2[df2$`#chrom1` %in% valid_chroms & df2$chrom2 %in% valid_chroms, ]

  df2 <- df2 %>%
    mutate(
      is_head1  = str_detect(source1, "\\d\\|"),
      is_head2  = str_detect(source2, "\\d\\|"),
      strand1   = ifelse(is_head1, "+", "-"),
      strand2   = ifelse(is_head2, "+", "-"),
      is_interchromosomal = (`#chrom1` != chrom2)
    ) %>%
    select(-is_head1, -is_head2, -is_interchromosomal)

  df3 <- df3[df3$seqnames %in% valid_chroms, ]

  # --- Build ShatterSeek input objects ---
  SV_data <- SVs(
    chrom1  = as.character(df2$`#chrom1`),
    pos1    = as.numeric(df2$start1),
    chrom2  = as.character(df2$chrom2),
    pos2    = as.numeric(df2$end2),
    SVtype  = as.character(df2$svtype),
    strand1 = as.character(df2$strand1),
    strand2 = as.character(df2$strand2)
  )

  CN_data <- CNVsegs(
    chrom    = as.character(df3$seqnames),
    start    = df3$start,
    end      = df3$end,
    total_cn = df3$TCN
  )

  # --- Run ShatterSeek ---
  chromothripsis <- shatterseek(SV.sample = SV_data, seg.sample = CN_data, genome = "hg19")
  summary        <- as.data.frame(chromothripsis@chromSummary)

  # --- Confidence classification ---
  val <- which(
    summary$clusterSize_including_TRA >= 6 &
    summary$max_number_oscillating_CN_segments_2_states >= 7 &
    summary$pval_fragment_joins < 0.05 &
    (summary$chr_breakpoint_enrichment < 0.05 | summary$pval_exp_chr < 0.05)
  )
  val <- paste0("High confidence  Chr:", val)

  val2 <- which(
    summary$clusterSize_including_TRA >= 3 &
    summary$number_TRA >= 4 &
    summary$pval_fragment_joins < 0.05 &
    summary$max_number_oscillating_CN_segments_2_states >= 7
  )
  val2 <- paste0("High confidence  Chr:", val2)

  val3 <- which(
    summary$clusterSize_including_TRA >= 6 &
    summary$max_number_oscillating_CN_segments_2_states >= 4 &
    summary$max_number_oscillating_CN_segments_2_states < 7 &
    summary$pval_fragment_joins < 0.05 &
    (summary$chr_breakpoint_enrichment <= 0.05 | summary$pval_exp_chr < 0.05)
  )
  val3 <- paste0("Low-confidence  Chr:", val3)

  max_length <- max(length(val), length(val2), length(val3))
  val  <- c(val,  rep(NA, max_length - length(val)))
  val2 <- c(val2, rep(NA, max_length - length(val2)))
  val3 <- c(val3, rep(NA, max_length - length(val3)))

  data <- data.frame(confidence1 = val, confidence2 = val2, confidence3 = val3)

  write.table(data,    file.path(output_dir, paste0(basename_i, ".chr_txt")))
  write.table(summary, file.path(output_dir, paste0(basename_i, ".summary_txt")),
              quote = FALSE, sep = "\t")
  qsave(chromothripsis, file.path(output_dir, paste0(basename_i, ".qs")))

  print(paste("-----------------------------------------------------------"))
}


##plotting the data 
data_dir <- "/path/to/output/shatterseek_results/"
file_list <- list.files(path = data_dir, pattern = "\\.chr_txt$", full.names = TRUE)



df_list <- lapply(file_list, read.table)
names(df_list) <- file_list
df <- df_list %>%
  rbindlist(idcol = 'sample_id', fill = T)

df_tidy <- df %>%
  # clean sample name from file path
  mutate(sample_name = basename(tools::file_path_sans_ext(sample_id))) %>%
  
  # pivot confidence columns to long format
  pivot_longer(
    cols      = starts_with("confidence"),
    names_to  = "conf_col",
    values_to = "label"
  ) %>%
  filter(!is.na(label)) %>%
  
  # parse label into confidence level and chromosome
  mutate(
    confidence = str_extract(label, "(?i)high confidence|low-confidence"),
    chr        = str_extract(label, "(?<=Chr:)\\S+"),
    chr        = ifelse(chr == "" | is.na(chr), NA_character_, paste0("chr", chr))
  ) %>%
  filter(!is.na(chr)) %>%
  select(conf_col, label, chr, confidence, sample_name)



# 4. Data Aggregation & Categorical Sorting for Plotting
plot_data <- df_tidy %>%
  group_by(chr, confidence) %>%
  tally(name = "Freq")

plot_data$chr <- gsub("chr", "",plot_data$chr)
chrom_order <- sort(as.numeric(unique(plot_data$chr)))
plot_data$chr <- factor(plot_data$chr, levels = as.character(chrom_order))

# 5. Generate and Save Stacked Bar Chart to PDF
pdf_output_path <- "chromosome_frequency_plot.pdf"
pdf(file = pdf_output_path, width = 10, height = 6)

p <- ggplot(plot_data, aes(x = chr, y = Freq, fill = confidence)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("High confidence" = "dodgerblue4", "Low-confidence" = "skyblue")) +
  theme_classic() +
  labs(
    title = "Chromosome Frequency Across Samples by Confidence Level",
    x = "Chromosome Number",
    y = "Total Occurrences",
    fill = "Confidence Status"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "top"
  )

print(p)
dev.off()


#-------------  Chromothripsis plot  ---------- #


plot_chromothripsis_grid <- function(
  chromothripsis,
  chromosomes,
  sample_name,
  genome = "hg19",
  heights = c(0.2, 0.4, 0.4, 0.4),
  ncol = NULL
) {
  # Build one arranged plot per chromosome
  chr_plots <- lapply(chromosomes, function(chr) {
    plots <- plot_chromothripsis(
      ShatterSeek_output = chromothripsis,
      chr               = chr,
      sample_name       = sample_name,
      genome            = genome
    )
    arrangeGrob(
      plots[[1]], plots[[2]], plots[[3]], plots[[4]],
      nrow    = 4,
      ncol    = 1,
      heights = heights,
      top     = chr   # optional: label each panel with its chromosome
    )
  })

  # Determine grid layout
  n <- length(chr_plots)
  ncol_final <- if (!is.null(ncol)) ncol else n

  # Draw the final combined plot
  do.call(plot_grid, c(chr_plots, list(ncol = ncol_final)))
}


setwd("/path/to/output/hrd/")
chromothripsis <- qread("ShatterSeek/analysis_output /svs_H021-5MKAFL_tumor-buffy_coat02_filtered_dedup_somatic.qs")
# Your original two-chromosome case
p <- plot_chromothripsis_grid(
  chromothripsis = chromothripsis,
  chromosomes    = c("chr4", "chr9"),
  sample_name    = "5MKAFL",
  genome         = "hg19"
)


pdf_output_path <- "ShatterSeek/shatter_5MKAFL.pdf"
pdf(file = pdf_output_path, width = 12, height = 6)
print(p)
dev.off()
