

log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "Figure6.log")
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


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


required_packages <- c("ggplot2", "dplyr", "tidyr", "ggpubr", "patchwork")


readr_loaded <- FALSE
try({

  suppressPackageStartupMessages(library(readr))
  readr_loaded <- TRUE
  log_message("readr package loaded successfully
", print_to_console = TRUE)
}, silent = TRUE)


if (!readr_loaded) {
  log_message("Failed to load readr, will use base read.csv
", print_to_console = TRUE)
}


for (pkg in required_packages) {

  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
    log_message(paste0("Installing package: ", pkg, "\n"), print_to_console = TRUE)
    install.packages(pkg)
  } else {
    log_message(paste0("Loading installed package: ", pkg, "\n"))
  }

  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


dplyr_version <- as.character(packageVersion("dplyr"))
log_message(paste0("dplyr version used: ", dplyr_version, "\n"), print_to_console = TRUE)


log_message("Loading previously generated GSVA scoring data...
", print_to_console = TRUE)


results_paths <- c(
  "results",
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


log_message("Loading TP53 and CTNNB1 mutation data...
", print_to_console = TRUE)


mutation_paths <- c(
  "data",
  "../data",
  file.path(getwd(), "data"),
  file.path(dirname(getwd()), "data")
)


existing_mutation_path <- NULL
for (path in mutation_paths) {
  if (dir.exists(path) &&
      file.exists(file.path(path, "alterations_across_samples_TP53.csv")) &&
      file.exists(file.path(path, "alterations_across_samples_CTNNB1.csv"))) {
    existing_mutation_path <- path
    log_message(paste0("Found mutation data file directory: ", existing_mutation_path, "\n"), print_to_console = TRUE)
    break
  }
}


if (is.null(existing_mutation_path)) {
  stop("Error: mutation data file not found. Please ensure it exists in the 'data' directory.")
}


tp53_mutations <- read.csv(file.path(existing_mutation_path, "alterations_across_samples_TP53.csv"))
ctnnb1_mutations <- read.csv(file.path(existing_mutation_path, "alterations_across_samples_CTNNB1.csv"))


log_message(paste0("TP53 mutation data column names: ", paste(colnames(tp53_mutations), collapse = ", "), "\n"))
log_message(paste0("CTNNB1 mutation data column names: ", paste(colnames(ctnnb1_mutations), collapse = ", "), "\n"))


log_message("Processing mutation data, creating mutation status metadata...
", print_to_console = TRUE)


if (!"Sample.ID" %in% colnames(tp53_mutations) || !"Sample.ID" %in% colnames(ctnnb1_mutations)) {

  sample_id_cols_tp53 <- grep("Sample|ID|Patient", colnames(tp53_mutations), ignore.case = TRUE)
  sample_id_cols_ctnnb1 <- grep("Sample|ID|Patient", colnames(ctnnb1_mutations), ignore.case = TRUE)

  if (length(sample_id_cols_tp53) > 0 && length(sample_id_cols_ctnnb1) > 0) {

    colnames(tp53_mutations)[sample_id_cols_tp53[1]] <- "Sample.ID"
    colnames(ctnnb1_mutations)[sample_id_cols_ctnnb1[1]] <- "Sample.ID"
    log_message("Renamed sample ID column to 'Sample.ID'
", print_to_console = TRUE)
  } else {
    stop("Error: sample ID column not found in mutation data.")
  }
}


altered_col_tp53 <- which(colnames(tp53_mutations) == "Altered")
altered_col_ctnnb1 <- which(colnames(ctnnb1_mutations) == "Altered")

if (length(altered_col_tp53) == 0 || length(altered_col_ctnnb1) == 0) {
  stop("Error: Altered column not found in mutation data.")
}

log_message(paste0("Using TP53 Altered column: ", colnames(tp53_mutations)[altered_col_tp53], "\n"))
log_message(paste0("Using CTNNB1 Altered column: ", colnames(ctnnb1_mutations)[altered_col_ctnnb1], "\n"))


mutation_data <- data.frame(
  PatientID = tp53_mutations$Sample.ID,
  TP53_Status = ifelse(tp53_mutations[[altered_col_tp53]] == 1, "Mutant", "WT"),
  stringsAsFactors = FALSE
)


ctnnb1_status <- data.frame(
  PatientID = ctnnb1_mutations$Sample.ID,
  CTNNB1_Status = ifelse(ctnnb1_mutations[[altered_col_ctnnb1]] == 1, "Mutant", "WT"),
  stringsAsFactors = FALSE
)


mutation_data <- merge(mutation_data, ctnnb1_status, by = "PatientID", all = TRUE)


mutation_data$TP53_Status[is.na(mutation_data$TP53_Status)] <- "Unknown"
mutation_data$CTNNB1_Status[is.na(mutation_data$CTNNB1_Status)] <- "Unknown"


log_message("TP53 mutation status statistics:
")
tp53_table <- table(mutation_data$TP53_Status)
log_message(paste(capture.output(print(tp53_table)), collapse = "\n"))
log_message("\n")

log_message("CTNNB1 mutation status statistics:
")
ctnnb1_table <- table(mutation_data$CTNNB1_Status)
log_message(paste(capture.output(print(ctnnb1_table)), collapse = "\n"))
log_message("\n")


log_message("Preparing long-format data...
", print_to_console = TRUE)


prepare_long_data_mutation <- function(gsva_scores, mutation_data) {

  gsva_long <- as.data.frame(t(gsva_scores))
  gsva_long$PatientID <- rownames(gsva_long)


  gsva_long$PatientID_short <- substr(gsva_long$PatientID, 1, 12)
  mutation_data$PatientID_short <- substr(mutation_data$PatientID, 1, 12)


  gsva_long_tidy <- tryCatch({

    gsva_long %>%
      pivot_longer(cols = !c(PatientID, PatientID_short),
                   names_to = "AgentState",
                   values_to = "Score")
  }, error = function(e) {

    log_message(paste0("pivot_longer failed, trying gather. Error: ", e$message, "\n"))

    tryCatch({
      gsva_long %>%
        gather(key = "AgentState", value = "Score",
               -PatientID, -PatientID_short)
    }, error = function(e2) {

      log_message(paste0("gather also failed, using base R. Error: ", e2$message, "\n"))


      cols_to_keep <- c("PatientID", "PatientID_short")
      cols_to_reshape <- setdiff(colnames(gsva_long), cols_to_keep)


      melted_data <- data.frame(
        PatientID = rep(gsva_long$PatientID, length(cols_to_reshape)),
        PatientID_short = rep(gsva_long$PatientID_short, length(cols_to_reshape)),
        AgentState = rep(cols_to_reshape, each = nrow(gsva_long)),
        Score = unlist(lapply(cols_to_reshape, function(col) gsva_long[[col]]))
      )

      return(melted_data)
    })
  })


  merged_data <- tryCatch({
    merge(gsva_long_tidy,
          mutation_data,
          by = "PatientID_short",
          all.x = TRUE)
  }, error = function(e) {
    log_message(paste0("merge failed, trying alternative. Error: ", e$message, "\n"))


    result <- merge(gsva_long_tidy, mutation_data,
                  by.x = "PatientID_short", by.y = "PatientID_short",
                  all.x = TRUE)
    return(result)
  })


  merged_data$TP53_Status[is.na(merged_data$TP53_Status)] <- "Unknown"
  merged_data$CTNNB1_Status[is.na(merged_data$CTNNB1_Status)] <- "Unknown"


  merged_data <- tryCatch({

    merged_data %>%
      filter(TP53_Status != "Unknown" & CTNNB1_Status != "Unknown")
  }, error = function(e) {
    log_message(paste0("dplyr filter failed, using base R. Error: ", e$message, "\n"))


    merged_data[merged_data$TP53_Status != "Unknown" &
                  merged_data$CTNNB1_Status != "Unknown", ]
  })


  merged_data$TP53_Status <- factor(merged_data$TP53_Status, levels = c("WT", "Mutant"))
  merged_data$CTNNB1_Status <- factor(merged_data$CTNNB1_Status, levels = c("WT", "Mutant"))

  return(merged_data)
}


long_data_tcga_mutation <- prepare_long_data_mutation(gsva_scores_tcga, mutation_data)


log_message(paste0("Long-format data dimensions: ", paste(dim(long_data_tcga_mutation), collapse = " x "), "\n"))
log_message(paste0("Samples in long-format data: ", length(unique(long_data_tcga_mutation$PatientID)), "\n"))


log_message("Selecting agent states to plot...
", print_to_console = TRUE)


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


all_states <- signature_data$Module_Display_Name
log_message(paste0("Total loaded: ", length(all_states), " signatures
"), print_to_console = TRUE)


tp53_related_states <- all_states


ctnnb1_related_states <- all_states


states_to_plot_mutation <- unique(c(tp53_related_states, ctnnb1_related_states))


states_to_plot_mutation <- intersect(states_to_plot_mutation, unique(long_data_tcga_mutation$AgentState))
if (length(states_to_plot_mutation) == 0) {
  stop("Error: selected agent states not in GSVA scoring data. Check state names.")
}
log_message(paste0("Will plot these agent states: ", paste(states_to_plot_mutation, collapse = ", "), "\n"), print_to_console = TRUE)


plot_data_mutation <- tryCatch({
  long_data_tcga_mutation %>%
    filter(AgentState %in% states_to_plot_mutation) %>%
    mutate(AgentState = factor(AgentState, levels = states_to_plot_mutation))
}, error = function(e) {
  log_message(paste0("dplyr filter/mutate failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


  result <- long_data_tcga_mutation[long_data_tcga_mutation$AgentState %in% states_to_plot_mutation, ]
  result$AgentState <- factor(result$AgentState, levels = states_to_plot_mutation)
  return(result)
})


log_message("Plotting TP53-related state figures...
", print_to_console = TRUE)


create_tp53_plot <- function(data, state_name) {

  state_data <- data %>% filter(AgentState == state_name)


  if (nrow(state_data) == 0) {
    log_message(paste0("Warning: ", state_name, " has no data
"), print_to_console = TRUE)

    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste("No data for", state_name)) +
      theme_minimal()
    return(p)
  }


  if (length(unique(state_data$TP53_Status)) < 2) {
    log_message(paste0("Warning: ", state_name, " has only one TP53 mutation group
"), print_to_console = TRUE)

    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste("Only one TP53 status group for", state_name)) +
      theme_minimal()
    return(p)
  }


  p_val <- tryCatch({
    wilcox.test(Score ~ TP53_Status, data = state_data)$p.value
  }, error = function(e) {
    log_message(paste0("Warning: Cannot compute", state_name, " p-value: ", e$message, "\n"), print_to_console = TRUE)
    return(NA)
  })


  p <- ggplot(state_data, aes(x = TP53_Status, y = Score, fill = TP53_Status)) +
    geom_violin(trim = FALSE, alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_manual(values = c("WT" = "lightblue", "Mutant" = "darkorange")) +
    labs(
      title = state_name,
      x = "TP53 Status",
      y = "Activity Score"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5, size = 9))


  if (!is.na(p_val)) {
    p <- p + stat_compare_means(
      aes(group = TP53_Status),
      method = "wilcox.test",
      label = "p.format",
      size = 3,
      label.y.npc = 0.9
    )
  }

  return(p)
}


log_message("Plotting CTNNB1-related state figures...
", print_to_console = TRUE)


create_ctnnb1_plot <- function(data, state_name) {

  state_data <- data %>% filter(AgentState == state_name)


  if (nrow(state_data) == 0) {
    log_message(paste0("Warning: ", state_name, " has no data
"), print_to_console = TRUE)

    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste("No data for", state_name)) +
      theme_minimal()
    return(p)
  }


  if (length(unique(state_data$CTNNB1_Status)) < 2) {
    log_message(paste0("Warning: ", state_name, " has only one CTNNB1 mutation group
"), print_to_console = TRUE)

    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste("Only one CTNNB1 status group for", state_name)) +
      theme_minimal()
    return(p)
  }


  p_val <- tryCatch({
    wilcox.test(Score ~ CTNNB1_Status, data = state_data)$p.value
  }, error = function(e) {
    log_message(paste0("Warning: Cannot compute", state_name, " p-value: ", e$message, "\n"), print_to_console = TRUE)
    return(NA)
  })


  p <- ggplot(state_data, aes(x = CTNNB1_Status, y = Score, fill = CTNNB1_Status)) +
    geom_violin(trim = FALSE, alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_manual(values = c("WT" = "lightgreen", "Mutant" = "purple")) +
    labs(
      title = state_name,
      x = "CTNNB1 Status",
      y = "Activity Score"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5, size = 9))


  if (!is.na(p_val)) {
    p <- p + stat_compare_means(
      aes(group = CTNNB1_Status),
      method = "wilcox.test",
      label = "p.format",
      size = 3,
      label.y.npc = 0.9
    )
  }

  return(p)
}


log_message("Saving data and figures...
", print_to_console = TRUE)


output_dir <- existing_results_path


log_message("Saving TP53 mutation status comparison data...
", print_to_console = TRUE)


tp53_means <- tryCatch({

  plot_data_mutation %>%
    group_by(AgentState, TP53_Status) %>%
    summarise(Mean_Score = mean(Score, na.rm = TRUE),
              SD_Score = sd(Score, na.rm = TRUE),
              N = n(),
              .groups = "drop")
}, error = function(e) {

  log_message(paste0("First summarise method failed, trying alternative. Error: ", e$message, "\n"), print_to_console = TRUE)
  tryCatch({
    plot_data_mutation %>%
      group_by(AgentState, TP53_Status) %>%
      summarise(Mean_Score = mean(Score, na.rm = TRUE),
                SD_Score = sd(Score, na.rm = TRUE),
                N = n()) %>%
      ungroup()
  }, error = function(e2) {

    log_message(paste0("Second summarise method also failed, using base R. Error: ", e2$message, "\n"), print_to_console = TRUE)
    result <- aggregate(Score ~ AgentState + TP53_Status, data = plot_data_mutation,
                      FUN = function(x) c(Mean_Score = mean(x, na.rm = TRUE),
                                        SD_Score = sd(x, na.rm = TRUE),
                                        N = length(x)))

    result_df <- do.call(data.frame, result)
    return(result_df)
  })
})


tp53_pvals <- data.frame(
  AgentState = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (state in unique(plot_data_mutation$AgentState)) {
  state_data <- plot_data_mutation[plot_data_mutation$AgentState == state, ]
  if (length(unique(state_data$TP53_Status)) >= 2) {
    p_val <- tryCatch({
      wilcox.test(Score ~ TP53_Status, data = state_data)$p.value
    }, error = function(e) {
      NA
    })
    tp53_pvals <- rbind(tp53_pvals, data.frame(AgentState = state, P_Value = p_val))
  }
}


tp53_means_wide <- tryCatch({
  tp53_means %>%
    pivot_wider(id_cols = AgentState,
                names_from = TP53_Status,
                values_from = c(Mean_Score, SD_Score, N))
}, error = function(e) {
  log_message(paste0("pivot_wider failed, trying alternative. Error: ", e$message, "\n"), print_to_console = TRUE)


  tryCatch({

    mean_data <- tp53_means %>%
      select(AgentState, TP53_Status, Mean_Score) %>%
      spread(key = TP53_Status, value = Mean_Score, sep = "_")

    sd_data <- tp53_means %>%
      select(AgentState, TP53_Status, SD_Score) %>%
      spread(key = TP53_Status, value = SD_Score, sep = "_")

    n_data <- tp53_means %>%
      select(AgentState, TP53_Status, N) %>%
      spread(key = TP53_Status, value = N, sep = "_")


    result <- left_join(mean_data, sd_data, by = "AgentState")
    result <- left_join(result, n_data, by = "AgentState")

    return(result)
  }, error = function(e2) {
    log_message(paste0("spread also failed, using base R. Error: ", e2$message, "\n"), print_to_console = TRUE)


    tp53_long <- tp53_means
    result <- reshape(tp53_long,
                    idvar = "AgentState",
                    timevar = "TP53_Status",
                    direction = "wide")


    colnames(result) <- gsub("Mean_Score\\.", "Mean_Score_", colnames(result))
    colnames(result) <- gsub("SD_Score\\.", "SD_Score_", colnames(result))
    colnames(result) <- gsub("N\\.", "N_", colnames(result))

    return(result)
  })
})


tp53_results <- left_join(tp53_means_wide, tp53_pvals, by = "AgentState")


tp53_csv <- file.path(output_dir, "Figure6_TP53_Mutation_Comparison.csv")
write.csv(tp53_results, tp53_csv, row.names = FALSE)
log_message(paste0("TP53 mutation comparison data table saved: ", tp53_csv, "\n"), print_to_console = TRUE)


log_message("Saving CTNNB1 mutation status comparison data...
", print_to_console = TRUE)


ctnnb1_means <- tryCatch({

  plot_data_mutation %>%
    group_by(AgentState, CTNNB1_Status) %>%
    summarise(Mean_Score = mean(Score, na.rm = TRUE),
              SD_Score = sd(Score, na.rm = TRUE),
              N = n(),
              .groups = "drop")
}, error = function(e) {

  log_message(paste0("First summarise method failed, trying alternative. Error: ", e$message, "\n"), print_to_console = TRUE)
  tryCatch({
    plot_data_mutation %>%
      group_by(AgentState, CTNNB1_Status) %>%
      summarise(Mean_Score = mean(Score, na.rm = TRUE),
                SD_Score = sd(Score, na.rm = TRUE),
                N = n()) %>%
      ungroup()
  }, error = function(e2) {

    log_message(paste0("Second summarise method also failed, using base R. Error: ", e2$message, "\n"), print_to_console = TRUE)
    result <- aggregate(Score ~ AgentState + CTNNB1_Status, data = plot_data_mutation,
                      FUN = function(x) c(Mean_Score = mean(x, na.rm = TRUE),
                                        SD_Score = sd(x, na.rm = TRUE),
                                        N = length(x)))

    result_df <- do.call(data.frame, result)
    return(result_df)
  })
})


ctnnb1_pvals <- data.frame(
  AgentState = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (state in unique(plot_data_mutation$AgentState)) {
  state_data <- plot_data_mutation[plot_data_mutation$AgentState == state, ]
  if (length(unique(state_data$CTNNB1_Status)) >= 2) {
    p_val <- tryCatch({
      wilcox.test(Score ~ CTNNB1_Status, data = state_data)$p.value
    }, error = function(e) {
      NA
    })
    ctnnb1_pvals <- rbind(ctnnb1_pvals, data.frame(AgentState = state, P_Value = p_val))
  }
}


ctnnb1_means_wide <- tryCatch({
  ctnnb1_means %>%
    pivot_wider(id_cols = AgentState,
                names_from = CTNNB1_Status,
                values_from = c(Mean_Score, SD_Score, N))
}, error = function(e) {
  log_message(paste0("pivot_wider failed, trying alternative. Error: ", e$message, "\n"), print_to_console = TRUE)


  tryCatch({

    mean_data <- ctnnb1_means %>%
      select(AgentState, CTNNB1_Status, Mean_Score) %>%
      spread(key = CTNNB1_Status, value = Mean_Score, sep = "_")

    sd_data <- ctnnb1_means %>%
      select(AgentState, CTNNB1_Status, SD_Score) %>%
      spread(key = CTNNB1_Status, value = SD_Score, sep = "_")

    n_data <- ctnnb1_means %>%
      select(AgentState, CTNNB1_Status, N) %>%
      spread(key = CTNNB1_Status, value = N, sep = "_")


    result <- left_join(mean_data, sd_data, by = "AgentState")
    result <- left_join(result, n_data, by = "AgentState")

    return(result)
  }, error = function(e2) {
    log_message(paste0("spread also failed, using base R. Error: ", e2$message, "\n"), print_to_console = TRUE)


    ctnnb1_long <- ctnnb1_means
    result <- reshape(ctnnb1_long,
                    idvar = "AgentState",
                    timevar = "CTNNB1_Status",
                    direction = "wide")


    colnames(result) <- gsub("Mean_Score\\.", "Mean_Score_", colnames(result))
    colnames(result) <- gsub("SD_Score\\.", "SD_Score_", colnames(result))
    colnames(result) <- gsub("N\\.", "N_", colnames(result))

    return(result)
  })
})


ctnnb1_results <- left_join(ctnnb1_means_wide, ctnnb1_pvals, by = "AgentState")


ctnnb1_csv <- file.path(output_dir, "Figure6_CTNNB1_Mutation_Comparison.csv")
write.csv(ctnnb1_results, ctnnb1_csv, row.names = FALSE)
log_message(paste0("CTNNB1 mutation comparison data table saved: ", ctnnb1_csv, "\n"), print_to_console = TRUE)


log_message("Outputting Type_II_IFN_Reponse pathway results...
", print_to_console = TRUE)


if (!"Type_II_IFN_Reponse" %in% unique(plot_data_mutation$AgentState)) {
  log_message("Warning: Type_II_IFN_Reponse pathway data not found
", print_to_console = TRUE)
} else {

  ifn_tp53_plot <- create_tp53_plot(plot_data_mutation, "Type_II_IFN_Reponse")


  ifn_ctnnb1_plot <- create_ctnnb1_plot(plot_data_mutation, "Type_II_IFN_Reponse")


  ifn_tp53_pdf <- file.path(output_dir, "Figure6E_TP53_Type_II_IFN_Reponse.pdf")
  ifn_ctnnb1_pdf <- file.path(output_dir, "Figure6E_CTNNB1_Type_II_IFN_Reponse.pdf")


  ggsave(ifn_tp53_pdf, ifn_tp53_plot, width = 4, height = 3, dpi = 300)
  ggsave(ifn_ctnnb1_pdf, ifn_ctnnb1_plot, width = 4, height = 3, dpi = 300)

  log_message(paste0("Type_II_IFN_Reponse TP53 figure saved: ", ifn_tp53_pdf, "\n"), print_to_console = TRUE)
  log_message(paste0("Type_II_IFN_Reponse CTNNB1 figure saved: ", ifn_ctnnb1_pdf, "\n"), print_to_console = TRUE)
}

log_message("Type_II_IFN_Reponse pathway analysis complete!
", print_to_console = TRUE)
