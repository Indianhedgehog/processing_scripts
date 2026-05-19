################################################################
# Adapted from ASCAT algorithm 
# plots overall as well as per chromosome CNV profile.
# plots the coverage as well as the copy number profile of all the chromosomes
# This workflow is adapted to be used with the output of the HMF pipeline, but can be adapted to be used with other CNV callers as well.
###################################################################################

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Provide complete PID along with tumor and control info \n")
}


pid <- args[1]
#pid <- "H021-1N79LK7_metastasis02_blood"

##Overall profile 
ascat.plotAdjustedAscatProfile=function(data_frame, REF, y_limit=5, plot_unrounded=FALSE, png_prefix="", purity = "", ploidy = "") {
  required_columns <- c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "sample")
  stopifnot(all(required_columns %in% colnames(data_frame)))
  
  data_frame <- data_frame[data_frame$chromosome != "Y", ]
  
  if (plot_unrounded) {
    SEGMENTS <- data_frame[, c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "sample")]
    colourA = "#943CC3" # purple
    colourB = "#60AF36" # green
  } else {
    SEGMENTS <- data_frame
    SEGMENTS$majorAlleleCopyNumber <- SEGMENTS$majorAlleleCopyNumber - 0.1
    SEGMENTS$minorAlleleCopyNumber <- SEGMENTS$minorAlleleCopyNumber + 0.1
    colourA = "#E03546" # red
    colourB = "#3557E0" # blue
  }
  SEGMENTS$majorAlleleCopyNumber=ifelse(SEGMENTS$majorAlleleCopyNumber>y_limit, y_limit+0.1, SEGMENTS$majorAlleleCopyNumber)
  SEGMENTS$minorAlleleCopyNumber=ifelse(SEGMENTS$minorAlleleCopyNumber>y_limit, y_limit+0.1, SEGMENTS$minorAlleleCopyNumber)
  
  if (REF=="hg19") {
    REF=data.frame(chrom=c(1:22, "X"),
                   start=rep(1, 23),
                   end=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                         135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                         59128983, 63025520, 48129895, 51304566, 155270560))
  } else if (REF=="hg38") {
    REF=data.frame(chrom=c(1:22, "X"),
                   start=rep(1, 23),
                   end=c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                         133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                         58617616, 64444167, 46709983, 50818468, 156040895))
  } else {
    stopifnot(is.data.frame(REF))
    stopifnot(identical(colnames(REF), c("chrom", "start", "end")))
  }
  
  SEGMENTS$chromosome=gsub("^chr", "", SEGMENTS$chromosome)
  stopifnot(all(unique(SEGMENTS$chromosome) %in% REF$chrom))
  REF$size=REF$end-REF$start+1
  REF$middle=0
  for (i in 1:nrow(REF)) {
    if (i==1) {
      REF$middle[i]=REF$size[i]/2
    } else {
      REF$middle[i]=sum(as.numeric(REF$size[1:(i-1)]))+REF$size[i]/2
    }
  }; rm(i)
  REF$cumul=cumsum(as.numeric(REF$size))
  REF$add=cumsum(as.numeric(c(0, REF$size[1:(nrow(REF)-1)])))
  
  SEGMENTS$startpos_adjusted=SEGMENTS$start
  SEGMENTS$endpos_adjusted=SEGMENTS$end
  for (CHR in unique(REF$chrom)) {
    INDEX=which(SEGMENTS$chromosome==CHR)
    if (length(INDEX)>0) {
      SEGMENTS$startpos_adjusted[INDEX]=SEGMENTS$startpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
      SEGMENTS$endpos_adjusted[INDEX]=SEGMENTS$endpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
    }
    rm(INDEX)
  }; rm(CHR)
  
  for (SAMPLE in sort(unique(SEGMENTS$sample))) {
    SEGS=SEGMENTS[which(SEGMENTS$sample==SAMPLE), ]
    if (nrow(SEGS)==0) warning(paste0("No segments for sample: ", SAMPLE))
    #maintitle = paste("Ploidy: ", sprintf("%1.2f", ASCAT_output_object$ploidy[SAMPLE]), ", purity: ", sprintf("%2.0f", ASCAT_output_object$purity[SAMPLE]*100), "%, goodness of fit: ", sprintf("%2.1f", ASCAT_output_object$goodnessOfFit[SAMPLE]), "%", ifelse(isTRUE(ASCAT_output_object$nonaberrantarrays[SAMPLE]), ", non-aberrant", ""), sep="")
    maintitle = paste("Ploidy: ", sprintf("%1.2f", ploidy), ", purity: ", sprintf("%2.0f", purity*100),"%")
    png(filename = paste0(png_prefix, SAMPLE, ".adjusted", ifelse(plot_unrounded, "rawprofile", "HMFprofile"), ".png"), width = 2000, height = (y_limit*100), res = 200)
    par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
    ticks=seq(0, y_limit, 1)
    plot(c(1, REF$cumul[nrow(REF)]), c(0, y_limit), type = "n", xaxt = "n", yaxt="n", main = maintitle, xlab = "", ylab = "")
    axis(side = 2, at = ticks)
    abline(h=ticks, col="lightgrey", lty=1)
    rect(SEGS$startpos_adjusted, (SEGS$majorAlleleCopyNumber-0.07), SEGS$endpos_adjusted, (SEGS$majorAlleleCopyNumber+0.07), col=ifelse(SEGS$majorAlleleCopyNumber>=y_limit, adjustcolor(colourA, red.f=0.75, green.f=0.75, blue.f=0.75), colourA), border=ifelse(SEGS$majorAlleleCopyNumber>=y_limit, adjustcolor(colourA, red.f=0.75, green.f=0.75, blue.f=0.75), colourA))
    rect(SEGS$startpos_adjusted, (SEGS$minorAlleleCopyNumber-0.07), SEGS$endpos_adjusted, (SEGS$minorAlleleCopyNumber+0.07), col=ifelse(SEGS$minorAlleleCopyNumber>=y_limit, adjustcolor(colourB, red.f=0.75, green.f=0.75, blue.f=0.75), colourB), border=ifelse(SEGS$minorAlleleCopyNumber>=y_limit, adjustcolor(colourB, red.f=0.75, green.f=0.75, blue.f=0.75), colourB))
    abline(v=c(1, REF$cumul), lty=1, col="lightgrey")
    text(REF$middle, y_limit, REF$chrom, pos = 1, cex = 2)
    dev.off()
    rm(SEGS, ticks, maintitle)
  }; rm(SAMPLE)
}

##By chromosome 
ascat.plotByChromosome <- function(data_frame, REF, y_limit=5, plot_unrounded=FALSE, png_prefix="", purity = "", ploidy = "") {
  required_columns <- c("chromosome", "start", "end", "majorAlleleCopyNumber", "minorAlleleCopyNumber", "sample")
  stopifnot(all(required_columns %in% colnames(data_frame)))
  
  # Remove Y chromosome
  data_frame <- data_frame[data_frame$chromosome != "Y", ]
  
  if (plot_unrounded) {
    colourA = "#943CC3" # purple
    colourB = "#60AF36" # green
  } else {
    data_frame$majorAlleleCopyNumber <- data_frame$majorAlleleCopyNumber - 0.1
    data_frame$minorAlleleCopyNumber <- data_frame$minorAlleleCopyNumber + 0.1
    colourA = "#E03546" # red
    colourB = "#3557E0" # blue
  }
  
  # Cap values at y_limit
  data_frame$majorAlleleCopyNumber <- pmin(data_frame$majorAlleleCopyNumber, y_limit + 0.1)
  data_frame$minorAlleleCopyNumber <- pmin(data_frame$minorAlleleCopyNumber, y_limit + 0.1)
  
  # Define chromosome reference
  if (REF == "hg19") {
    REF <- data.frame(chrom = c(1:22, "X"),
                      start = rep(1, 23),
                      end = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                              135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                              59128983, 63025520, 48129895, 51304566, 155270560))
  } else if (REF == "hg38") {
    REF <- data.frame(chrom = c(1:22, "X"),
                      start = rep(1, 23),
                      end = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                              133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                              58617616, 64444167, 46709983, 50818468, 156040895))
  } else {
    stopifnot(is.data.frame(REF))
    stopifnot(identical(colnames(REF), c("chrom", "start", "end")))
  }
  
  data_frame$chromosome <- gsub("^chromosome", "", data_frame$chromosome)
  stopifnot(all(unique(data_frame$chromosome) %in% REF$chrom))
  
  # Loop through each sample and chromosome
  for (SAMPLE in unique(data_frame$sample)) {
    for (CHR in unique(data_frame$chromosome)) {
      SEGS <- data_frame[data_frame$sample == SAMPLE & data_frame$chromosome == CHR, ]
      if (nrow(SEGS) == 0) next  # Skip if no segments
      
      chrom_size <- REF$end[REF$chrom == CHR]
      
      maintitle <- paste("Chromosome:", CHR, "Ploidy:", sprintf("%1.2f", ploidy), "Purity:", sprintf("%2.0f", purity * 100), "%")
      png(filename = paste0(png_prefix, SAMPLE, "_chr", CHR, "_profile.png"), width = 1500, height = 800, res = 200)
      
      par(mar = c(4, 5, 5, 2), cex = 0.8, cex.main = 1.5, cex.axis = 1.2)
      plot(c(1, chrom_size), c(0, y_limit), type = "n", xaxt = "n", yaxt = "n", main = maintitle, xlab = "Genomic Position", ylab = "Copy Number")
      axis(side = 1, at = seq(0, chrom_size, length.out = 5), labels = round(seq(0, chrom_size / 1e6, length.out = 5), 1))
      axis(side = 2, at = seq(0, y_limit, 1))
      abline(h = seq(0, y_limit, 1), col = "lightgrey", lty = 1)
      
      # Plot segments
      rect(SEGS$start, SEGS$majorAlleleCopyNumber - 0.07, SEGS$end, SEGS$majorAlleleCopyNumber + 0.07, col = colourA, border = colourA)
      rect(SEGS$start, SEGS$minorAlleleCopyNumber - 0.07, SEGS$end, SEGS$minorAlleleCopyNumber + 0.07, col = colourB, border = colourB)
      
      dev.off()
    }
  }
}


cnv_path <- file.path(pid, "purple", paste0(pid, "_T.purple.cnv.somatic.tsv"))
cnv_file <- readr::read_tsv(cnv_path)
cnv_file$sample <- pid

pp_path <- file.path(pid, "purple", paste0(pid, "_T.purple.purity.tsv"))
pp_file <- readr::read_tsv(pp_path)

##create the dir 
plots_dir <- file.path(pid, "plots/cnv_profile")

if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
}

## Output files
png_prefix <- file.path(plots_dir, "output_")

ascat.plotAdjustedAscatProfile(cnv_file, REF = "hg19", y_limit = 5, plot_unrounded = F, png_prefix = png_prefix, purity = pp_file$purity , ploidy = pp_file$ploidy)

ascat.plotByChromosome(cnv_file, REF = "hg19", y_limit = 5, plot_unrounded = F, png_prefix = png_prefix, purity = pp_file$purity, ploidy = pp_file$ploidy)



###################################################################
# Log2 ratio's from the tumor/normal files
# Depends on cobalt ratio's
#################################################################

# Input parameters
cobalt_file <- file.path(pid, "cobalt", paste0(pid, "_T.cobalt.ratio.tsv.gz"))
sex <- pp_file$gender 

plots_dir <- file.path(pid, "plots/coverage")

if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
}


# Read COBALT data
coverage_data <- read_tsv(cobalt_file)
# Expected columns: chromosome, position, tumorRatio (or similar)

# Define chromosomes
#chr_count <- ifelse(sex == "MALE", 24, 23)
#chr_names <- if(sex == "MALE") c(1:22, 'X', 'Y') else c(1:22, 'X')


# Prepare all three data types
coverage_all <- coverage_data %>%
  mutate(
    tumor_log2 = log2(ifelse(tumorGCRatio > 0, tumorGCRatio, NA)),
    normal_log2 = log2(ifelse(referenceGCRatio > 0, referenceGCRatio, NA)),
    ratio_log2 = log2(ifelse(tumorGCRatio > 0 & referenceGCRatio > 0, 
                             tumorGCRatio / referenceGCRatio, NA))
  ) %>%
  pivot_longer(cols = c(normal_log2, tumor_log2, ratio_log2),
               names_to = "plot_type",
               values_to = "log2_value") %>%
  filter(is.finite(log2_value)) %>%
  mutate(plot_type = factor(plot_type,
                            levels = c("normal_log2", "tumor_log2", "ratio_log2"),
                            labels = c("Normal", "Tumor", "Tumor/Normal Ratio")))

# Get chromosomes
chromosomes <- unique(coverage_all$chromosome)
chr_order <- if(sex == "MALE") c(1:22, "X", "Y") else c(1:22, "X")
chromosomes <- chromosomes[order(match(chromosomes, chr_order))]

# Loop through each chromosome
for (chr in chromosomes) {
  
  chr_data <- coverage_all %>% filter(chromosome == chr)
  
  # Create three-panel plot
  p <- ggplot(chr_data, aes(x = position/1e6, y = log2_value)) +
    geom_point(size = 0.3, alpha = 0.5, 
               aes(color = plot_type)) +
    geom_hline(yintercept = 0, col = "red", lwd = 0.6) +
    geom_hline(yintercept = c(-2, -1, 1, 2), col = "gray", lty = "dashed", lwd = 0.3) +
    facet_grid(plot_type ~ ., scales = "fixed") +
    scale_color_manual(values = c("Normal" = "blue3", 
                                   "Tumor" = "red3", 
                                   "Tumor/Normal Ratio" = "steelblue")) +
    ylim(-4, 4) +
    labs(x = "Position (MB)", 
         y = "log2 Value",
         title = paste0(sample_id, " - Chromosome ", chr, ", sex = ", sex)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "gray90", color = NA),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10)
    )
  
  # Save plot
  filename <- file.path(output_dir, paste0("chr", chr, "_coverage_3panel_10kb.pdf"))
  ggsave(filename, p, width = 12, height = 10, dpi = 300)
  
#  cat(paste0("Saved: ", filename, "\n"))
}



# Aggregate to 10kb windows
coverage_10kb <- coverage_data %>%
  mutate(
    # Create 10kb bins (position divided by 10000, then rounded down)
    window_10kb = floor(position / 10000) * 10000
  ) %>%
  group_by(chromosome, window_10kb) %>%
  summarise(
    # Average the values across the 10 1kb windows
    tumorGCRatio = mean(tumorGCRatio, na.rm = TRUE),
    referenceGCRatio = mean(referenceGCRatio, na.rm = TRUE),
    # Keep the window start position
    position = first(window_10kb),
    .groups = 'drop'
  )

# Now use coverage_10kb for plotting
coverage_all <- coverage_10kb %>%
  mutate(
    tumor_log2 = log2(ifelse(tumorGCRatio > 0, tumorGCRatio, NA)),
    normal_log2 = log2(ifelse(referenceGCRatio > 0, referenceGCRatio, NA)),
    ratio_log2 = log2(ifelse(tumorGCRatio > 0 & referenceGCRatio > 0, 
                             tumorGCRatio / referenceGCRatio, NA))
  ) %>%
  pivot_longer(cols = c(normal_log2, tumor_log2, ratio_log2),
               names_to = "plot_type",
               values_to = "log2_value") %>%
  filter(is.finite(log2_value)) %>%
  mutate(plot_type = factor(plot_type,
                            levels = c("normal_log2", "tumor_log2", "ratio_log2"),
                            labels = c("Normal", "Tumor", "Tumor/Normal Ratio")))

# Continue with your plotting code...
