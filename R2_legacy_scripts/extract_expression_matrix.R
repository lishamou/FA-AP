

suppressPackageStartupMessages(library(dplyr))


log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "extract_expression_matrix.log")
log_con <- file(log_file, "w")
log_message <- function(message, print_to_console = FALSE) {
  cat(message, file = log_con, append = TRUE)
  if (print_to_console) cat(message)
}
on.exit({
  close(log_con)
  log_message(paste0("\nLog saved to: ", normalizePath(log_file), "\n"), print_to_console = TRUE)
})


log_message("Reading TCGA-LIHC gene expression data...
", print_to_console = TRUE)


data_paths <- c(
  "data/TCGA-LIHC_expression_data.txt",
  "../data/TCGA-LIHC_expression_data.txt",
  file.path(getwd(), "data/TCGA-LIHC_expression_data.txt"),
  file.path(dirname(getwd()), "data/TCGA-LIHC_expression_data.txt")
)


existing_path <- NULL
for (path in data_paths) {
  if (file.exists(path)) {
    existing_path <- path
    log_message(paste0("Found data file: ", existing_path, "\n"), print_to_console = TRUE)
    break
  }
}
if (is.null(existing_path)) {
  stop("Error: TCGA-LIHC_expression_data.txt not found. Please ensure the file exists in the 'data' directory.")
}


log_message("Reading expression data and checking duplicate row names...
")
temp_data <- read.table(existing_path, header = TRUE, sep = "\t", check.names = FALSE)


first_col_name <- colnames(temp_data)[1]
first_col_values <- temp_data[[first_col_name]]
if (any(duplicated(first_col_values))) {
  log_message("Warning: duplicate row identifiers found, will use unique row names
")
  dup_rows <- first_col_values[duplicated(first_col_values)]
  log_message(paste0("Duplicate row identifiers: ", paste(dup_rows, collapse = ", "), "\n"))
  unique_rownames <- make.unique(as.character(first_col_values))
  rownames(temp_data) <- unique_rownames
  expr_data <- temp_data[, -1, drop = FALSE]
} else {
  log_message("No duplicate row identifiers found
")
  rownames(temp_data) <- first_col_values
  expr_data <- temp_data[, -1, drop = FALSE]
}


if (ncol(expr_data) > nrow(expr_data)) {
  expr_data <- t(expr_data)
}


log_message("Converting data frame to matrix...
")
expr_data <- as.matrix(expr_data)


if (max(expr_data, na.rm = TRUE) > 50) {
  expr_matrix <- log2(expr_data + 1)
  log_message("Applied log2(TPM+1) conversion to expression data
")
} else {
  expr_matrix <- expr_data
  log_message("Expression data already log scale, no conversion needed
")
}


log_message("Starting TCGA sample filtering...
")

select_tumor_rna_tcga_samples <- function(sample_barcodes) {
  if (length(sample_barcodes) == 0) return(character(0))
  barcode_df <- data.frame(full_barcode = sample_barcodes, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      sample_type_code = substr(full_barcode, 14, 15),
      patient_sample_type = substr(full_barcode, 1, 15),
      vial = substr(full_barcode, 16, 16),
      analyte = substr(full_barcode, 20, 20),
      plate = substr(full_barcode, 22, 25),
      patient_id = substr(full_barcode, 1, 12)
    )
  tumor_barcode_df <- barcode_df %>%
    dplyr::filter(as.numeric(sample_type_code) >= 1 & as.numeric(sample_type_code) <= 9)
  if (nrow(tumor_barcode_df) == 0) return(character(0))
  selected_df <- tumor_barcode_df %>%
    dplyr::group_by(patient_sample_type) %>%
    dplyr::filter(if (any(vial == 'A')) vial == 'A' else TRUE) %>%
    dplyr::mutate(analyte_priority = dplyr::case_when(analyte == 'R' ~ 1, analyte == 'T' ~ 2, TRUE ~ 99)) %>%
    dplyr::filter(analyte_priority == min(analyte_priority)) %>%
    dplyr::filter(analyte_priority != 99) %>%
    dplyr::filter(n() > 0) %>%
    dplyr::arrange(desc(plate)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  duplicated_patients <- selected_df$patient_id[duplicated(selected_df$patient_id)]
  if (length(duplicated_patients) > 0) {
    log_message(sprintf("Found %d duplicate patient IDs: %s\n", length(duplicated_patients), paste(duplicated_patients, collapse = ", ")))
    selected_df <- selected_df %>%
      dplyr::group_by(patient_id) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    log_message("Processed duplicate samples, one sample per patient
")
  }
  return(selected_df$full_barcode)
}

sample_names <- colnames(expr_matrix)
tcga_pattern <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[0-9]{2}"
is_tcga_sample <- grepl(tcga_pattern, sample_names)
if (any(is_tcga_sample)) {
  log_message("Detected TCGA samples, starting filtering...
")
  tcga_samples <- sample_names[is_tcga_sample]
  selected_samples <- select_tumor_rna_tcga_samples(tcga_samples)
  if (length(selected_samples) > 0) {
    expr_matrix <- expr_matrix[, selected_samples, drop = FALSE]
    log_message(sprintf("Sample filtering done. Original count: %d, , filtered count: %d\n", length(tcga_samples), length(selected_samples)))
  } else {
    stop("Error: no samples after filtering")
  }
} else {
  log_message("Warning: TCGA-format sample names not detected, skipping filtering step
")
}


dir.create("results", showWarnings = FALSE)
expr_matrix_file <- file.path("results", "TCGA-LIHC_expression_matrix_processed.rds")
saveRDS(expr_matrix, expr_matrix_file)
log_message(paste0("Processed expression matrix saved: ", expr_matrix_file, "\n"), print_to_console = TRUE)
