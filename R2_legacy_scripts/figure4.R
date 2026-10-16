

log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "Figure4.log")
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


required_packages <- c("GSVA", "ggplot2", "dplyr", "tidyr", "readxl", "openxlsx",
                      "reshape2", "ggpubr", "patchwork")


for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  log_message(paste0("Loading package: ", pkg, "\n"))
}


log_message("Loading previously generated GSVA scoring data...
", print_to_console = TRUE)


results_paths <- c(
  "results",
  "./results",
  "../results",
  file.path(getwd(), "results"),
  file.path(dirname(getwd()), "results")
)


existing_results_path <- NULL
for (path in results_paths) {
  if (dir.exists(path) &&
      file.exists(file.path(path, "gsva_scores.rds")) &&
      file.exists(file.path(path, "patient_metadata.rds"))) {
    existing_results_path <- path
    log_message(paste0("Found data file directory: ", existing_results_path, "\n"), print_to_console = TRUE)
    break
  }
}


if (is.null(existing_results_path)) {
  stop("Error: previously generated data file not found. Please run figure2.R first.")
}


gsva_scores_tcga <- readRDS(file.path(existing_results_path, "gsva_scores.rds"))
patient_metadata <- readRDS(file.path(existing_results_path, "patient_metadata.rds"))


log_message("Reading risk score data...
", print_to_console = TRUE)


risk_score_paths <- c(
  "data/The risk score of each patient.xlsx",
  "../data/The risk score of each patient.xlsx",
  file.path(getwd(), "data/The risk score of each patient.xlsx"),
  file.path(dirname(getwd()), "data/The risk score of each patient.xlsx")
)


existing_risk_score_path <- NULL
for (path in risk_score_paths) {
  if (file.exists(path)) {
    existing_risk_score_path <- path
    log_message(paste0("Found risk score data file: ", existing_risk_score_path, "\n"), print_to_console = TRUE)
    break
  }
}


if (is.null(existing_risk_score_path)) {
  stop("Error: 'The risk score of each patient.xlsx' not found. Please ensure the file exists in the 'data' directory.")
}


risk_score_data <- read.xlsx(existing_risk_score_path)


log_message(paste0("Risk score data column names: ", paste(colnames(risk_score_data), collapse = ", "), "\n"))


log_message("Assigning risk groups to TCGA samples...
", print_to_console = TRUE)


sample_id_col <- NULL
risk_score_col <- NULL


for (col in colnames(risk_score_data)) {
  if (any(grepl("TCGA-", risk_score_data[[col]]))) {
    sample_id_col <- col
    log_message(paste0("Using sample ID column from risk score data: ", sample_id_col, "\n"))
    break
  }
}


possible_score_cols <- c("risk_score", "Risk_Score", "RS", "score", "Score")
for (col in colnames(risk_score_data)) {
  if (col %in% possible_score_cols || grepl("score|Score|risk|Risk", col)) {
    risk_score_col <- col
    log_message(paste0("Using risk score column from risk score data: ", risk_score_col, "\n"))
    break
  }
}


if (is.null(sample_id_col)) {
  sample_id_col <- colnames(risk_score_data)[1]
  log_message(paste0("Sample ID column not found, using first column: ", sample_id_col, "\n"))
}
if (is.null(risk_score_col)) {
  risk_score_col <- colnames(risk_score_data)[2]
  log_message(paste0("Risk score column not found, using second column: ", risk_score_col, "\n"))
}


risk_data <- data.frame(
  SampleID = risk_score_data[[sample_id_col]],
  RiskScore = as.numeric(risk_score_data[[risk_score_col]])
)


risk_data$SampleID_short <- substr(risk_data$SampleID, 1, 12)


risk_median <- median(risk_data$RiskScore, na.rm = TRUE)
log_message(paste0("Risk score median: ", risk_median, "\n"))


risk_data$RiskGroup <- ifelse(risk_data$RiskScore > risk_median, "High", "Low")
log_message(paste0("High-risk group sample count: ", sum(risk_data$RiskGroup == "High", na.rm = TRUE), "\n"), print_to_console = TRUE)
log_message(paste0("Low-risk group sample count: ", sum(risk_data$RiskGroup == "Low", na.rm = TRUE), "\n"), print_to_console = TRUE)


log_message("Matching risk group info with GSVA scoring data...
")


gsva_samples_short <- substr(colnames(gsva_scores_tcga), 1, 12)


match_indices <- match(gsva_samples_short, risk_data$SampleID_short)


risk_groups <- risk_data$RiskGroup[match_indices]
names(risk_groups) <- colnames(gsva_scores_tcga)


risk_groups[is.na(risk_groups)] <- "Unknown"


high_risk_samples <- colnames(gsva_scores_tcga)[risk_groups == "High"]
low_risk_samples <- colnames(gsva_scores_tcga)[risk_groups == "Low"]

log_message(paste0("High-risk samples in GSVA data: ", length(high_risk_samples), "\n"), print_to_console = TRUE)
log_message(paste0("Low-risk samples in GSVA data: ", length(low_risk_samples), "\n"), print_to_console = TRUE)


log_message("Preparing plot data...
")


log_message("Reading all signatures from signature.csv...
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


prepare_long_data <- function(gsva_scores, high_samples, low_samples, cohort_name) {

  high_risk_scores <- gsva_scores[, high_samples, drop = FALSE]
  low_risk_scores <- gsva_scores[, low_samples, drop = FALSE]


  high_risk_long <- reshape2::melt(high_risk_scores, varnames = c("AgentState", "SampleID"),
                                  value.name = "Score")
  high_risk_long$RiskGroup <- "High"
  high_risk_long$Cohort <- cohort_name

  low_risk_long <- reshape2::melt(low_risk_scores, varnames = c("AgentState", "SampleID"),
                                 value.name = "Score")
  low_risk_long$RiskGroup <- "Low"
  low_risk_long$Cohort <- cohort_name


  combined_long <- rbind(high_risk_long, low_risk_long)


  combined_long$RiskGroup <- factor(combined_long$RiskGroup, levels = c("Low", "High"))


  return(combined_long)
}


tcga_long_data <- prepare_long_data(gsva_scores_tcga, high_risk_samples, low_risk_samples, "TCGA-LIHC")


log_message("Saving data...
")


output_dir <- existing_results_path


log_message("Saving inter-group risk data...
")


risk_group_means <- tcga_long_data %>%
  group_by(AgentState, RiskGroup) %>%
  summarize(Mean_Score = mean(Score, na.rm = TRUE),
            SD_Score = sd(Score, na.rm = TRUE),
            N = n(),
            .groups = "drop")


risk_group_pvals <- data.frame(
  AgentState = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (state in unique(tcga_long_data$AgentState)) {
  state_data <- tcga_long_data[tcga_long_data$AgentState == state, ]
  if (length(unique(state_data$RiskGroup)) >= 2) {
    p_val <- tryCatch({
      wilcox.test(Score ~ RiskGroup, data = state_data)$p.value
    }, error = function(e) {
      NA
    })
    risk_group_pvals <- rbind(risk_group_pvals, data.frame(AgentState = state, P_Value = p_val))
  }
}


risk_group_means_wide <- risk_group_means %>%
  pivot_wider(id_cols = AgentState,
              names_from = RiskGroup,
              values_from = c(Mean_Score, SD_Score, N))


risk_group_results <- left_join(risk_group_means_wide, risk_group_pvals, by = "AgentState")


risk_group_csv <- file.path(output_dir, "Figure4_Risk_Group_Comparison.csv")
write.csv(risk_group_results, risk_group_csv, row.names = FALSE)
log_message(paste0("Inter-group risk data table saved: ", risk_group_csv, "\n"))


long_data_csv <- file.path(output_dir, "Figure4_Risk_Group_Long_Format.csv")
write.csv(tcga_long_data, long_data_csv, row.names = FALSE)
log_message(paste0("Risk group long-format data table saved: ", long_data_csv, "\n"))
