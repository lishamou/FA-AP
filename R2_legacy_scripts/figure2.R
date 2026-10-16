

log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "figure2.log")
log_con <- file(log_file, "w")


log_message <- function(message, print_to_console = FALSE) {
  cat(message, file = log_con, append = TRUE)
  if (print_to_console) {
    cat(message)
  }
}


on.exit({
  close(log_con)
  log_message(paste0("\nLog saved to: ", normalizePath(log_file), "\n"), print_to_console = TRUE)
})


log_message("Creating output directory...
", print_to_console = TRUE)

results_paths <- c(
  "results",
  "../results",
  file.path(getwd(), "results"),
  file.path(dirname(getwd()), "results")
)


output_dir <- NULL
for (path in results_paths) {

  if (dir.exists(path)) {
    output_dir <- path
    log_message(paste0("Using existing output directory: ", output_dir, "\n"))
    break
  }

  tryCatch({
    dir.create(path, recursive = TRUE)
    output_dir <- path
    log_message(paste0("Creating new output directory: ", output_dir, "\n"))
    break
  }, error = function(e) {

    log_message(paste0("Cannot create directory: ", path, "\n"))
  })
}


if (is.null(output_dir)) {
  stop("Error: cannot create output directory. Please ensure you have write permission.")
}


required_packages <- c("GSVA", "pheatmap", "RColorBrewer", "dplyr",
                      "readxl", "openxlsx")


for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  log_message(paste0("Loading package: ", pkg, "\n"))
}


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
  if (length(sample_barcodes) == 0) {
    return(character(0))
  }


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
    dplyr::filter(

      as.numeric(sample_type_code) >= 1 & as.numeric(sample_type_code) <= 9
    )


  if (nrow(tumor_barcode_df) == 0) {
    return(character(0))
  }


  selected_df <- tumor_barcode_df %>%

    dplyr::group_by(patient_sample_type) %>%

    dplyr::filter(
      if (any(vial == 'A')) vial == 'A' else TRUE
    ) %>%

    dplyr::mutate(
      analyte_priority = dplyr::case_when(
        analyte == 'R' ~ 1,
        analyte == 'T' ~ 2,
        TRUE ~ 99
      )
    ) %>%
    dplyr::filter(analyte_priority == min(analyte_priority)) %>%
    dplyr::filter(analyte_priority != 99) %>%

    dplyr::filter(n() > 0) %>%
    dplyr::arrange(desc(plate)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()


  duplicated_patients <- selected_df$patient_id[duplicated(selected_df$patient_id)]

  if (length(duplicated_patients) > 0) {

    log_message(sprintf("Found %d duplicate patient IDs: %s\n",
                       length(duplicated_patients),
                       paste(duplicated_patients, collapse = ", ")))


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
    log_message(sprintf("Sample filtering done. Original count: %d, , filtered count: %d\n",
                       length(tcga_samples), length(selected_samples)))
  } else {
    stop("Error: no samples after filtering")
  }
} else {
  log_message("Warning: TCGA-format sample names not detected, skipping filtering step
")
}


log_message("Reading TCGA-LIHC clinical data...
", print_to_console = TRUE)


clinical_paths <- c(
  "data/TCGA-LIHC_clinical_data.xlsx",
  "../data/TCGA-LIHC_clinical_data.xlsx",
  file.path(getwd(), "data/TCGA-LIHC_clinical_data.xlsx"),
  file.path(dirname(getwd()), "data/TCGA-LIHC_clinical_data.xlsx")
)


existing_clinical_path <- NULL
for (path in clinical_paths) {
  if (file.exists(path)) {
    existing_clinical_path <- path
    log_message(paste0("Found clinical data file: ", existing_clinical_path, "\n"), print_to_console = TRUE)
    break
  }
}


if (is.null(existing_clinical_path)) {
  stop("Error: TCGA-LIHC_clinical_data.xlsx not found. Please ensure it exists in the 'data' directory.")
}


clinical_data <- read.xlsx(existing_clinical_path)


log_message("Processing clinical and prognosis data...
")


log_message("Converting sample IDs in clinical data to short format...
")
clinical_data$patient_id_short <- substr(clinical_data$patient_id, 1, 12)


colnames_short <- substr(colnames(expr_matrix), 1, 12)


log_message("Matching clinical data with expression matrix...
")
matched_indices <- match(colnames_short, clinical_data$patient_id_short)


total_samples_before_match <- ncol(expr_matrix)
log_message(sprintf("\nTotal samples before matching: %d\n", total_samples_before_match))


patient_metadata <- data.frame(
  row.names = colnames(expr_matrix)
)


patient_metadata$PatientID_short <- substr(rownames(patient_metadata), 1, 12)


patient_metadata$SurvivalStatus <- clinical_data$fustat[matched_indices]


patient_metadata$OS <- clinical_data$futime[matched_indices]


patient_metadata$TP53_Status <- clinical_data$TP53_Status[matched_indices]


patient_metadata$CTNNB1_Status <- clinical_data$CTNNB1_Status[matched_indices]


patient_metadata$RiskGroup <- clinical_data$Risk_Group[matched_indices]


patient_metadata <- patient_metadata[!is.na(patient_metadata$SurvivalStatus) &
                                   !is.na(patient_metadata$TP53_Status) &
                                   !is.na(patient_metadata$CTNNB1_Status) &
                                   !is.na(patient_metadata$RiskGroup), ]


patient_metadata$RiskGroup <- ifelse(tolower(patient_metadata$RiskGroup) == "high", "High", "Low")


total_samples_after_match <- nrow(patient_metadata)
missing_samples <- total_samples_before_match - total_samples_after_match

log_message(sprintf("\nSample matching and data integrity statistics:
"))
log_message(sprintf("Samples before matching: %d\n", total_samples_before_match))
log_message(sprintf("Complete samples after matching: %d\n", total_samples_after_match))
log_message(sprintf("Samples excluded due to incomplete data: %d\n", missing_samples))
log_message(sprintf("Data integrity rate: %.2f%%\n", total_samples_after_match/total_samples_before_match*100))


log_message("\nClinical features of final analyzed samples:
")


clinical_stats <- list()

log_message("\n1. Survival status distribution:
")
survival_table <- table(patient_metadata$SurvivalStatus)
log_message(paste(capture.output(print(survival_table)), collapse = "\n"))
survival_rate <- survival_table["Alive"]/sum(survival_table)*100
log_message(sprintf("\nSurvival rate: %.2f%%\n", survival_rate))


clinical_stats[["Survival_Status"]] <- data.frame(
  Category = "Survival Status",
  Group = names(survival_table),
  Count = as.numeric(survival_table),
  Percentage = sprintf("%.2f%%", as.numeric(survival_table)/sum(survival_table)*100)
)

log_message("\n2. TP53 mutation status distribution:
")
tp53_table <- table(patient_metadata$TP53_Status)
log_message(paste(capture.output(print(tp53_table)), collapse = "\n"))
tp53_rate <- tp53_table["Mutant"]/sum(tp53_table)*100
log_message(sprintf("\nTP53 mutation rate: %.2f%%\n", tp53_rate))


clinical_stats[["TP53_Status"]] <- data.frame(
  Category = "TP53 Status",
  Group = names(tp53_table),
  Count = as.numeric(tp53_table),
  Percentage = sprintf("%.2f%%", as.numeric(tp53_table)/sum(tp53_table)*100)
)

log_message("\n3. CTNNB1 mutation status distribution:
")
ctnnb1_table <- table(patient_metadata$CTNNB1_Status)
log_message(paste(capture.output(print(ctnnb1_table)), collapse = "\n"))
ctnnb1_rate <- ctnnb1_table["Mutant"]/sum(ctnnb1_table)*100
log_message(sprintf("\nCTNNB1 mutation rate: %.2f%%\n", ctnnb1_rate))


clinical_stats[["CTNNB1_Status"]] <- data.frame(
  Category = "CTNNB1 Status",
  Group = names(ctnnb1_table),
  Count = as.numeric(ctnnb1_table),
  Percentage = sprintf("%.2f%%", as.numeric(ctnnb1_table)/sum(ctnnb1_table)*100)
)

log_message("\n4. Risk group distribution:
")
risk_table <- table(patient_metadata$RiskGroup)
log_message(paste(capture.output(print(risk_table)), collapse = "\n"))
high_risk_rate <- risk_table["High"]/sum(risk_table)*100
log_message(sprintf("\nHigh-risk proportion: %.2f%%\n", high_risk_rate))


clinical_stats[["Risk_Group"]] <- data.frame(
  Category = "Risk Group",
  Group = names(risk_table),
  Count = as.numeric(risk_table),
  Percentage = sprintf("%.2f%%", as.numeric(risk_table)/sum(risk_table)*100)
)


clinical_stats_df <- do.call(rbind, clinical_stats)
rownames(clinical_stats_df) <- NULL


sample_stats <- data.frame(
  Category = c("Sample Statistics", "Sample Statistics", "Sample Statistics", "Sample Statistics"),
  Group = c("Total Samples Before Match", "Complete Samples After Match", "Excluded Samples", "Data Completeness"),
  Count = c(total_samples_before_match, total_samples_after_match, missing_samples, NA),
  Percentage = c(NA, NA, NA, sprintf("%.2f%%", total_samples_after_match/total_samples_before_match*100))
)


final_stats_df <- rbind(sample_stats, clinical_stats_df)


clinical_stats_file <- file.path(output_dir, "Figure2A_Clinical_Features_Summary.csv")
write.csv(final_stats_df, clinical_stats_file, row.names = FALSE)
log_message(sprintf("\nClinical feature statistics saved to: %s\n", clinical_stats_file))


log_message("Reading actual signatures from signature.csv...
", print_to_console = TRUE)


signature_paths <- c(
  "data/signature.csv",
  "../data/signature.csv",
  file.path(getwd(), "data/signature.csv"),
  file.path(dirname(getwd()), "data/signature.csv")
)


existing_signature_path <- NULL
for (path in signature_paths) {
  if (file.exists(path)) {
    existing_signature_path <- path
    log_message(paste0("Found signature file: ", existing_signature_path, "\n"), print_to_console = TRUE)
    break
  }
}


if (is.null(existing_signature_path)) {
  stop("Error: signature.csv not found. Please ensure it exists in the 'data' directory.")
}


signature_data <- read.csv(existing_signature_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
log_message(paste0("Successfully read ", nrow(signature_data), " gene set signatures
"), print_to_console = TRUE)


agent_signatures_list <- list()


for (i in 1:nrow(signature_data)) {

  geneset_name <- signature_data$Module_Display_Name[i]
  genes_str <- signature_data$Gene_Names[i]


  genes <- unlist(strsplit(genes_str, ",\\s*"))


  genes <- trimws(genes)


  agent_signatures_list[[geneset_name]] <- genes

  log_message(paste0("Adding gene set: ", geneset_name, " (", length(genes), " genes)
"))
}

log_message(paste0("Total loaded: ", length(agent_signatures_list), " gene set signatures
"), print_to_console = TRUE)


agent_signatures_list <- lapply(agent_signatures_list, function(genes) {
  genes_present <- intersect(genes, rownames(expr_matrix))
  if (length(genes_present) < 3) {
    log_message(paste0("Warning: signature ", names(agent_signatures_list)[which(agent_signatures_list == genes)],
        " only has ", length(genes_present), " genes in expression matrix
"))
  }
  return(genes_present)
})


agent_signatures_list <- agent_signatures_list[sapply(agent_signatures_list, length) >= 3]


log_message("Computing GSVA activity scores...
", print_to_console = TRUE)
gsva_scores <- gsva(expr_matrix, agent_signatures_list, method = "gsva", verbose = TRUE)


log_message("Preparing heatmap annotations...
")


annotation_col <- patient_metadata[colnames(gsva_scores), c("RiskGroup", "SurvivalStatus", "TP53_Status", "CTNNB1_Status")]


ann_colors <- list(
  RiskGroup = c(High = "firebrick", Low = "steelblue"),
  SurvivalStatus = c(Alive = "forestgreen", Dead = "grey30"),
  TP53_Status = c(Mutant = "darkorange", WT = "lightblue"),
  CTNNB1_Status = c(Mutant = "purple", WT = "lightgreen")
)


log_message("Plotting agent activity spectrum heatmap...
", print_to_console = TRUE)


heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)


output_pdf <- file.path(output_dir, "Figure2A_Agent_Activity_Heatmap.pdf")
pdf(output_pdf, width = 12, height = 16)


heatmap_result <- pheatmap(
  gsva_scores,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 12,
  fontsize_col = 10,
  border_color = NA,
  main = "Agent Functional State Activity Profile in TCGA-LIHC",
  silent = TRUE
)


pheatmap(
  gsva_scores,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 12,
  fontsize_col = 10,
  border_color = NA,
  main = "Agent Functional State Activity Profile in TCGA-LIHC"
)

dev.off()
log_message(paste0("PDF file saved: ", output_pdf, "\n"), print_to_console = TRUE)


output_png <- file.path(output_dir, "Figure2A_Agent_Activity_Heatmap.png")
png(output_png, width = 1200, height = 1600, res = 120)
pheatmap(
  gsva_scores,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  fontsize_col = 6,
  border_color = NA,
  main = "Agent Functional State Activity Profile in TCGA-LIHC",
  cellwidth = NA,
  cellheight = NA,
  treeheight_row = 50,
  treeheight_col = 30,
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  legend = TRUE
)
dev.off()
log_message(paste0("PNG file saved: ", output_png, "\n"), print_to_console = TRUE)


log_message("Saving processed data...
")


gsva_scores_file <- file.path(output_dir, "gsva_scores.rds")
patient_metadata_file <- file.path(output_dir, "patient_metadata.rds")
saveRDS(gsva_scores, gsva_scores_file)
saveRDS(patient_metadata, patient_metadata_file)
log_message(paste0("RDS data file saved: ", gsva_scores_file, "\n"))
log_message(paste0("RDS data file saved: ", patient_metadata_file, "\n"))


gsva_scores_csv <- file.path(output_dir, "Figure2A_GSVA_Scores.csv")

gsva_scores_df <- as.data.frame(gsva_scores)
gsva_scores_df$GeneSet <- rownames(gsva_scores)

gsva_scores_df <- gsva_scores_df[, c(ncol(gsva_scores_df), 1:(ncol(gsva_scores_df)-1))]
write.csv(gsva_scores_df, gsva_scores_csv, row.names = FALSE)
log_message(paste0("Raw GSVA scoring data table saved: ", gsva_scores_csv, "\n"))


clustered_matrix <- gsva_scores

clustered_gsva_scores_csv <- file.path(output_dir, "Figure2A_GSVA_Scores_Clustered.csv")
clustered_gsva_scores_df <- as.data.frame(clustered_matrix)
clustered_gsva_scores_df$GeneSet <- rownames(clustered_matrix)
clustered_gsva_scores_df <- clustered_gsva_scores_df[, c(ncol(clustered_gsva_scores_df), 1:(ncol(clustered_gsva_scores_df)-1))]
write.csv(clustered_gsva_scores_df, clustered_gsva_scores_csv, row.names = FALSE)
log_message(paste0("Clustered GSVA scoring data table saved: ", clustered_gsva_scores_csv, "\n"))


clustering_info <- list(
  row_clustering = NULL,
  col_clustering = NULL
)
clustering_info_file <- file.path(output_dir, "Figure2A_Clustering_Info.rds")
saveRDS(clustering_info, clustering_info_file)
log_message(paste0("Clustering tree info saved: ", clustering_info_file, "\n"))


patient_metadata_csv <- file.path(output_dir, "Figure2A_Patient_Metadata.csv")

patient_metadata_df <- as.data.frame(patient_metadata)
patient_metadata_df$SampleID <- rownames(patient_metadata)

patient_metadata_df <- patient_metadata_df[, c(ncol(patient_metadata_df), 1:(ncol(patient_metadata_df)-1))]
write.csv(patient_metadata_df, patient_metadata_csv, row.names = FALSE)
log_message(paste0("Patient metadata table saved: ", patient_metadata_csv, "\n"))

log_message(paste0("Figure 2A plotting complete! Results saved in ", output_dir, " directory
"), print_to_console = TRUE)
