
# build_oncoanalyzer_metadata.R
# Usage:
#   Rscript build_oncoanalyzer_metadata.R \
#     --wgs_dir     /omics/odcf/.../whole_genome_sequencing/view-by-pid \
#     --rna_dir     /omics/odcf/.../rna_sequencing/view-by-pid \
#     --output      fastq_rna_dna_metadata.csv \
#     --samples_file chordoma_samples_raw.tsv \  # optional
#     --pids        H021-XXXXX,H021-YYYYY        # optional

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(stringr)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-s", "--samples_file"), type = "character", default = NULL,
              help = "Optional TSV with columns: ProjectID, PatientID, DNASeq.
                      If omitted, ALL subdirectories under wgs_dir/rna_dir are scanned."),
  make_option(c("-w", "--wgs_dir"), type = "character",
              help = "WGS view-by-pid base directory [required]"),
  make_option(c("-r", "--rna_dir"), type = "character", default = NULL,
              help = "RNA view-by-pid base directory [optional]"),
  make_option(c("-o", "--output"), type = "character",
              default = "fastq_rna_dna_metadata.csv",
              help = "Output CSV file path [default: fastq_rna_dna_metadata.csv]"),
  make_option(c("-p", "--pids"), type = "character", default = NULL,
              help = "Optional comma-separated PIDs to filter (e.g. H021-XXXXXX,H021-YYYYYY)")
)

# ── Interactive vs CLI ────────────────────────────────────────────────────────
if (interactive()) {
  message("[DEV] Running interactively — using hardcoded opt list")
  opt <- list(
    samples_file = NULL,   # set to NULL to scan all dirs, or provide a path
    wgs_dir      = "/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid",
    rna_dir      = "/omics/odcf/project/hipo/hipo_021/sequencing/rna_sequencing/view-by-pid",
    output       = "fastq_rna_dna_metadata_test.csv",
    pids         = NULL    # e.g. "H021-1E41KP,H021-2577B8"
  )
} else {
  opt <- parse_args(OptionParser(option_list = option_list))
}

# ── Validate required args ────────────────────────────────────────────────────
if (is.null(opt$wgs_dir)) stop("Missing required argument: --wgs_dir")
if (!dir.exists(opt$wgs_dir)) stop("wgs_dir not found: ", opt$wgs_dir)

# ── Utility: clean subject ID ─────────────────────────────────────────────────
clean_subject_id <- function(x) {
  sub("_buffy[-_]coat.*|_tumor.*|_metastasis.*|_blood.*|_control.*|_undefined-neoplasia.*|-rep.*",
      "", x, perl = TRUE)
}

# ── Utility: pair R1 + R2 fastqs ─────────────────────────────────────────────
pair_fastqs <- function(patient_ids, base_dir) {

  patient_ids <- patient_ids[!is.na(patient_ids) & nzchar(patient_ids)]
  if (length(patient_ids) == 0) {
    warning("pair_fastqs: no valid patient IDs provided.")
    return(tibble())
  }
  if (is.null(base_dir) || !nzchar(base_dir)) {
    stop("pair_fastqs: base_dir is NULL or empty.")
  }

  all_files <- unlist(lapply(patient_ids, function(pid) {
    path <- file.path(base_dir, pid)
    if (length(path) == 0 || is.na(path) || !nzchar(path)) {
      message("  [WARN] Skipping empty path for PID: ", pid)
      return(character(0))
    }
    if (!dir.exists(path)) {
      message("  [WARN] Directory not found, skipping: ", path)
      return(character(0))
    }
    list.files(path, pattern = "\\.fastq\\.gz$", recursive = TRUE, full.names = TRUE)
  })) %>% unique()

  if (length(all_files) == 0) {
    warning("pair_fastqs: no fastq.gz files found under: ", base_dir)
    return(tibble())
  }

  make_half <- function(read_tag) {
    files <- grep(read_tag, all_files, value = TRUE)
    if (length(files) == 0) return(tibble())
    tibble(
      fullpath   = files,
      pids       = sub(paste0("^", base_dir, "/?"), "", files) %>%
                   sub("/.*", "", .),
      sample_ids = basename(files) %>%
        str_remove("\\.fastq\\.gz$") %>%
        str_remove(paste0("_", read_tag))
    )
  }

  r1 <- make_half("R1")
  r2 <- make_half("R2")

  if (nrow(r1) == 0 || nrow(r2) == 0) {
    warning("pair_fastqs: could not find R1 or R2 files under: ", base_dir)
    return(tibble())
  }

  unmatched_r1 <- anti_join(r1, r2, by = "sample_ids")
  unmatched_r2 <- anti_join(r2, r1, by = "sample_ids")
  if (nrow(unmatched_r1) > 0)
    message("  [WARN] ", nrow(unmatched_r1), " R1 file(s) have no R2 match — skipped")
  if (nrow(unmatched_r2) > 0)
    message("  [WARN] ", nrow(unmatched_r2), " R2 file(s) have no R1 match — skipped")

  inner_join(r1, r2, by = "sample_ids", suffix = c("_R1", "_R2")) %>%
    transmute(
      sample_ids,
      pids     = pids_R1,
      filepath = paste0(fullpath_R1, ";", fullpath_R2)
    )
}

# ── Resolve patient IDs ───────────────────────────────────────────────────────
# Mode A: samples_file provided  → read PIDs from TSV
# Mode B: no samples_file        → discover PIDs from directory listing
message("\n[1/4] Resolving patient IDs...")

if (!is.null(opt$samples_file)) {
  message("      Mode: samples file — ", opt$samples_file)
  if (!file.exists(opt$samples_file)) stop("samples_file not found: ", opt$samples_file)

  files_df <- read_tsv(opt$samples_file, show_col_types = FALSE) %>%
    filter(DNASeq == "WGS") %>%
    mutate(PatientID = paste0(ProjectID, "-", PatientID)) %>%
    filter(!is.na(PatientID), nzchar(PatientID))

  patient_ids <- files_df$PatientID

} else {
  message("      Mode: auto-discover — scanning: ", opt$wgs_dir)
  patient_ids <- list.dirs(opt$wgs_dir, full.names = FALSE, recursive = FALSE)
  patient_ids <- patient_ids[nzchar(patient_ids)]
}

# Optional PID filter
if (!is.null(opt$pids)) {
  pid_filter  <- str_split(opt$pids, ",")[[1]] %>% str_trim()
  message("      Filtering to ", length(pid_filter), " PID(s): ",
          paste(pid_filter, collapse = ", "))
  patient_ids <- patient_ids[patient_ids %in% pid_filter]
  if (length(patient_ids) == 0) stop("No matching PIDs found after filtering.")
}

message("      Found ", length(patient_ids), " patient(s)")

# ── Patterns ──────────────────────────────────────────────────────────────────
TUMOR_PATTERN   <- "tumor|metastasis"
CONTROL_PATTERN <- "blood|buffy[-_]coat|control"
TUMOR_NAME_RE   <- "\\b(metastasis[[:alnum:]_\\-]*|tumor[[:alnum:]_\\-]*)\\b"
CONTROL_NAME_RE <- "\\b(blood[[:alnum:]_\\-]*|buffy[-_]coat[[:alnum:]_\\-]*|control[[:alnum:]_\\-]*)\\b"
RNA_NAME_RE     <- "\\b(metastasis[[:alnum:]_\\-]*|tumor[[:alnum:]_\\-]*|undefined[[:alnum:]_\\-]*)\\b"

# ── WGS DNA ───────────────────────────────────────────────────────────────────
message("\n[2/4] Processing WGS DNA fastqs from: ", opt$wgs_dir)

wgs_paired <- pair_fastqs(patient_ids, opt$wgs_dir)
if (nrow(wgs_paired) == 0) stop("No WGS fastq pairs found. Check --wgs_dir.")

tumor_wgs <- wgs_paired %>%
  filter(grepl(TUMOR_PATTERN, filepath, ignore.case = TRUE)) %>%
  mutate(
    names       = str_extract(filepath, TUMOR_NAME_RE),
    group_id    = paste0(pids, "_", names),
    sample_type = "tumor",
    pid_key     = pids
  )

ctrl_wgs <- wgs_paired %>%
  filter(grepl(CONTROL_PATTERN, filepath, ignore.case = TRUE)) %>%
  mutate(
    names       = str_extract(filepath, CONTROL_NAME_RE),
    sample_type = "normal",
    pid_key     = pids
  )

if (nrow(tumor_wgs) == 0) stop("No tumor WGS samples found (tumor|metastasis pattern).")
if (nrow(ctrl_wgs)  == 0) stop("No control WGS samples found (blood|buffy|control pattern).")

matched <- merge(tumor_wgs, ctrl_wgs, by = "pid_key", suffix = c("_t", "_c"))
if (nrow(matched) == 0) stop("No tumor-control pairs could be matched by patient ID.")

t_out <- matched %>% transmute(filepath = filepath_t, group_id, names = names_t, sample_type = "tumor")
c_out <- matched %>% transmute(filepath = filepath_c, group_id, names = names_c, sample_type = "normal")

wgs_bind <- bind_rows(t_out, c_out) %>%
  arrange(group_id) %>%
  distinct() %>%
  mutate(
    filetype      = "fastq",
    sequence_type = "dna",
    subject_id    = clean_subject_id(group_id),
    sample_id     = case_when(
      sample_type == "tumor"  ~ paste0(group_id, "_T"),
      sample_type == "normal" ~ paste0(group_id, "_N")
    ),
    info = paste0(
      "library_id:S1;",
      rep(c("lane:001", "lane:002"), length.out = n())
    )
  )

message("      WGS rows: ", nrow(wgs_bind),
        " (", sum(wgs_bind$sample_type == "tumor"),  " tumor, ",
              sum(wgs_bind$sample_type == "normal"), " normal)")

# ── RNA ───────────────────────────────────────────────────────────────────────
message("\n[3/4] Processing RNA fastqs...")

if (is.null(opt$rna_dir)) {
  message("      --rna_dir not provided — skipping RNA.")
  rna_df <- tibble()
} else if (!dir.exists(opt$rna_dir)) {
  message("      [WARN] rna_dir not found, skipping RNA: ", opt$rna_dir)
  rna_df <- tibble()
} else {
  message("      From: ", opt$rna_dir)
  rna_paired <- pair_fastqs(patient_ids, opt$rna_dir)

  if (nrow(rna_paired) == 0) {
    message("      [WARN] No RNA fastq pairs found — skipping.")
    rna_df <- tibble()
  } else {
    rna_df <- rna_paired %>%
      mutate(
        names         = str_extract(filepath, RNA_NAME_RE),
        group_id      = paste0(pids, "_", names),
        filetype      = "fastq",
        sequence_type = "rna",
        sample_type   = "tumor",
        subject_id    = clean_subject_id(group_id),
        sample_id     = paste0(subject_id, "_RNA"),
        info          = paste0(
          "library_id:S3;",
          rep(c("lane:001", "lane:002"), length.out = n())
        )
      ) %>%
      filter(subject_id %in% wgs_bind$subject_id)

    message("      RNA rows (matched to WGS subjects): ", nrow(rna_df))
  }
}

# ── Combine + deduplicate ─────────────────────────────────────────────────────
message("\n[4/4] Combining, deduplicating, writing output...")

cols <- c("filepath", "group_id", "filetype", "sequence_type",
          "sample_type", "subject_id", "sample_id", "info")

final <- bind_rows(
    wgs_bind[, cols],
    if (nrow(rna_df) > 0) rna_df[, intersect(cols, colnames(rna_df))] else tibble()
  ) %>%
  arrange(subject_id) %>%
  distinct(filepath, group_id, .keep_all = TRUE)

write.table(final, opt$output,
            quote = FALSE, row.names = FALSE, sep = ",")

message("[DONE] Written: ", opt$output,
        "\n       Rows : ", nrow(final),
        "\n       DNA  : ", sum(final$sequence_type == "dna"),
        "\n       RNA  : ", sum(final$sequence_type == "rna"))

invisible(final)


# =============================================================================
# OPTIONAL FUNCTION: write_chunk()
# =============================================================================
write_chunk <- function(metadata_csv, chunk_group_ids, output_csv) {
  if (!file.exists(metadata_csv)) stop("Metadata file not found: ", metadata_csv)

  chunk <- read_csv(metadata_csv, show_col_types = FALSE) %>%
    filter(group_id %in% chunk_group_ids) %>%
    mutate(
      info = paste0(
        "library_id:S", sprintf("%03d", row_number()),
        ";lane:",        sprintf("%03d", row_number())
      )
    )

  if (nrow(chunk) == 0) {
    warning("write_chunk: no rows matched the provided group_ids.")
    return(invisible(NULL))
  }

  write.table(chunk, output_csv, quote = FALSE, row.names = FALSE, sep = ",")
  message("[CHUNK] Written: ", output_csv, " (", nrow(chunk), " rows)")
  invisible(chunk)
}
'''

with open('/root/scripts/build_oncoanalyzer_metadata.R', 'w') as f:
    f.write(script)
print("Written successfully")