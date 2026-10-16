

log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "figure2B.log")
log_con <- file(log_file, "w")


Sys.setlocale("LC_ALL", "en_US.UTF-8")


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


required_packages <- c("GSVA", "pheatmap", "RColorBrewer", "dplyr",
                      "readxl", "openxlsx", "ggplot2", "reshape2")


for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  log_message(paste0("Loading package: ", pkg, "\n"))
}


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


gsva_scores <- readRDS(file.path(existing_results_path, "gsva_scores.rds"))
patient_metadata <- readRDS(file.path(existing_results_path, "patient_metadata.rds"))


log_message("Computing average activity in different groups...
", print_to_console = TRUE)


log_message("Extracting survival status from patient_metadata...
", print_to_console = TRUE)


clustered_data_path <- file.path(existing_results_path, "Figure2A_GSVA_Scores_Clustered.csv")
if (!file.exists(clustered_data_path)) {
  stop("Error: Figure2A_GSVA_Scores_Clustered.csv not found. Please run figure2.R first.")
}
clustered_data <- read.csv(clustered_data_path)
clustered_state_order <- clustered_data$GeneSet
log_message(paste0("Loaded clustered State order, total ", length(clustered_state_order), " states
"), print_to_console = TRUE)


if ("SurvivalStatus" %in% colnames(patient_metadata)) {
  log_message("Found SurvivalStatus column, processing survival status data...
", print_to_console = TRUE)


  log_message(paste0("Samples in GSVA data: ", ncol(gsva_scores), "\n"))
  log_message(paste0("Samples in patient_metadata: ", nrow(patient_metadata), "\n"))
  log_message(paste0("First 5 GSVA samples: ", paste(head(colnames(gsva_scores), 5), collapse = ", "), "\n"))
  log_message(paste0("First 5 patient_metadata samples: ", paste(head(rownames(patient_metadata), 5), collapse = ", "), "\n"))


  common_samples <- intersect(colnames(gsva_scores), rownames(patient_metadata))
  log_message(paste0("Found ", length(common_samples), " matched samples
"), print_to_console = TRUE)

  if (length(common_samples) > 0) {

    survival_status <- patient_metadata[common_samples, "SurvivalStatus"]
    alive_samples <- common_samples[survival_status == "Alive"]
    dead_samples <- common_samples[survival_status == "Dead"]

    log_message(paste0("Alive sample count: ", length(alive_samples), "\n"), print_to_console = TRUE)
    log_message(paste0("Dead sample count: ", length(dead_samples), "\n"), print_to_console = TRUE)

    if (length(alive_samples) > 0 && length(dead_samples) > 0) {

      alive_mean <- rowMeans(gsva_scores[, alive_samples, drop = FALSE], na.rm = TRUE)
      dead_mean <- rowMeans(gsva_scores[, dead_samples, drop = FALSE], na.rm = TRUE)


      survival_status_means <- data.frame(
        State = rownames(gsva_scores),
        Alive = alive_mean,
        Dead = dead_mean
      )


      survival_status_means$Difference <- survival_status_means$Dead - survival_status_means$Alive


      survival_status_means$State <- factor(survival_status_means$State,
                                          levels = clustered_state_order,
                                          ordered = TRUE)


      survival_status_means <- survival_status_means[order(survival_status_means$State), ]


      survival_status_means_long <- reshape2::melt(survival_status_means,
                                                  id.vars = "State",
                                                  measure.vars = c("Dead", "Alive"),
                                                  variable.name = "Survival_Status",
                                                  value.name = "Activity")


      survival_status_means_long <- survival_status_means_long[order(survival_status_means_long$State), ]


      log_message("Plotting bar plot...
")


      p1 <- ggplot(survival_status_means_long, aes(x = State, y = Activity, fill = Survival_Status)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(values = c("Dead" = "grey30", "Alive" = "forestgreen")) +
        theme_light() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          panel.border = element_rect(color = "grey90", fill = NA, linewidth = 0.5),
          plot.background = element_rect(fill = "white", color = NA)
        ) +
        labs(x = "Agent Functional State", y = "Mean GSVA Score",
             title = "Agent Functional State Activity by Survival Status (Horizontal)",
             fill = "Survival Status")


      p2 <- ggplot(survival_status_means_long, aes(x = Activity, y = factor(State, levels = rev(clustered_state_order)), fill = Survival_Status)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_fill_manual(values = c("Dead" = "grey30", "Alive" = "forestgreen")) +
        theme_light() +
        theme(
          axis.text.y = element_text(size = 8, color = "black", hjust = 0),
          axis.text.x = element_text(size = 8, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          legend.position = "top",
          plot.margin = margin(5, 20, 5, 5, "mm"),
          panel.border = element_rect(color = "grey90", fill = NA, linewidth = 0.5),
          plot.background = element_rect(fill = "white", color = NA),

          axis.text.y.right = element_text(size = 8, color = "black", hjust = 0),
          axis.title.y.right = element_text(size = 12, color = "black"),
          axis.text.y.left = element_blank(),
          axis.title.y.left = element_blank(),
          axis.ticks.y.left = element_blank()
        ) +
        scale_y_discrete(position = "right") +
        labs(y = "Agent Functional State", x = "Mean GSVA Score",
             title = "Agent Functional State Activity by Survival Status",
             fill = "Survival Status")


      output_pdf <- file.path(existing_results_path, "Figure2B_Survival_Status_Comparison_Horizontal.pdf")
      output_png <- file.path(existing_results_path, "Figure2B_Survival_Status_Comparison_Horizontal.png")
      ggsave(output_pdf, p1, width = 20, height = 6, bg = "white")
      ggsave(output_png, p1, width = 20, height = 6, dpi = 300, bg = "white")


      output_pdf_vertical <- file.path(existing_results_path, "Figure2B_Survival_Status_Comparison_Vertical.pdf")
      output_png_vertical <- file.path(existing_results_path, "Figure2B_Survival_Status_Comparison_Vertical.png")
      ggsave(output_pdf_vertical, p2, width = 8, height = 12, bg = "white")
      ggsave(output_png_vertical, p2, width = 8, height = 12, dpi = 300, bg = "white")


      survival_status_csv <- file.path(existing_results_path, "Figure2B_Survival_Status_Comparison.csv")
      write.csv(survival_status_means, survival_status_csv, row.names = FALSE)
      log_message(paste0("Survival status comparison data table saved: ", survival_status_csv, "\n"))


      survival_status_long_csv <- file.path(existing_results_path, "Figure2B_Survival_Status_Long_Format.csv")
      write.csv(survival_status_means_long, survival_status_long_csv, row.names = FALSE)
      log_message(paste0("Survival status long-format data table saved: ", survival_status_long_csv, "\n"))

      log_message("Survival status comparison figures (horizontal and vertical) generated
", print_to_console = TRUE)
    } else {
      log_message("Warning: insufficient alive or dead samples, skipping survival status comparison
", print_to_console = TRUE)
    }
  } else {
    log_message("Warning: no matched samples, skipping survival status comparison
", print_to_console = TRUE)
  }
} else {
  log_message("Warning: SurvivalStatus column not found in patient_metadata, skipping survival status comparison
", print_to_console = TRUE)
}


if (exists("p1") || exists("p2")) {
  log_message("Generating combined comparison figure...
", print_to_console = TRUE)


  suppressPackageStartupMessages(library(gridExtra))
  log_message("Loading package: gridExtra\n")


  plot_list <- list()
  if (exists("p1")) plot_list[["Survival Status (Horizontal)"]] <- p1
  if (exists("p2")) plot_list[["Survival Status (Vertical)"]] <- p2


  combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 1)


  output_pdf <- file.path(existing_results_path, "Figure2B_Combined_Comparisons.pdf")
  output_png <- file.path(existing_results_path, "Figure2B_Combined_Comparisons.png")
  ggsave(output_pdf, combined_plot, width = 20, height = 6 * length(plot_list))
  ggsave(output_png, combined_plot, width = 20, height = 6 * length(plot_list), dpi = 300)

  log_message("Combined comparison figure generated
", print_to_console = TRUE)
} else {
  log_message("Warning: no comparison figures generated, skipping combined figure
", print_to_console = TRUE)
}

log_message(paste0("Figure 2B plotting complete! Results saved in ", existing_results_path, " directory
"), print_to_console = TRUE)
