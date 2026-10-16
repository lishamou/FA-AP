

log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "figure2B_final_visualization.log")
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


log_message("Loading required R packages...
", print_to_console = TRUE)
required_packages <- c("ggplot2", "dplyr", "pheatmap", "RColorBrewer")


for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  log_message(paste0("Loading package: ", pkg, "\n"))
}


log_message("Searching for results directory...
", print_to_console = TRUE)
results_paths <- c(
  "results",
  "../results",
  file.path(getwd(), "results"),
  file.path(dirname(getwd()), "results")
)


existing_results_path <- NULL
for (path in results_paths) {
  if (dir.exists(path) && file.exists(file.path(path, "Figure2B_Survival_Status_Comparison.csv"))) {
    existing_results_path <- path
    log_message(paste0("Found data file directory: ", existing_results_path, "\n"), print_to_console = TRUE)
    break
  }
}


if (is.null(existing_results_path)) {
  stop("Error: Figure2B_Survival_Status_Comparison.csv not found. Please run figure2B.R first.")
}


survival_dir <- file.path(existing_results_path, "survival_analysis")
has_survival_analysis <- dir.exists(survival_dir) &&
                       file.exists(file.path(survival_dir, "Cox_Regression_Results.csv"))


log_message("Reading analysis data...
", print_to_console = TRUE)
comparison_file <- file.path(existing_results_path, "Figure2B_Survival_Status_Comparison.csv")
survival_data <- read.csv(comparison_file, stringsAsFactors = FALSE)


cox_results <- NULL
if (has_survival_analysis) {
  log_message("Reading Cox regression analysis results...
", print_to_console = TRUE)
  cox_results_file <- file.path(survival_dir, "Cox_Regression_Results.csv")
  cox_results <- read.csv(cox_results_file, stringsAsFactors = FALSE)
}


log_message("Creating and saving each figure...
", print_to_console = TRUE)


top_n_states <- 10
survival_data$AbsDiff <- abs(survival_data$Difference)
top_states <- survival_data[order(survival_data$AbsDiff, decreasing = TRUE), ][1:top_n_states, ]
top_states$Direction <- ifelse(top_states$Difference > 0, "Higher in Dead", "Higher in Alive")


top_states$State <- factor(top_states$State,
                          levels = rev(top_states$State[order(top_states$Difference)]))


panel1 <- ggplot(top_states, aes(x = State, y = Difference, fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Higher in Dead" = "red", "Higher in Alive" = "blue")) +
  coord_flip() +
  theme_light() +
  labs(
    title = "Functional State Activity Differences",
    y = "Difference (Dead - Alive)",
    x = "Functional State",
    fill = "Enrichment Group"
  ) +
  theme(
    plot.title = element_text(size = 9, face = "bold", color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 9, color = "black"),
    legend.text = element_text(size = 9, color = "black"),

    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.line.y.left = element_blank(),

    axis.text.y.right = element_text(size = 9, color = "black", hjust = 0),
    axis.title.y.right = element_text(size = 9, color = "black", angle = 90),
    axis.line.y.right = element_line(color = "black"),
    axis.ticks.y.right = element_line(color = "black"),

    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top = element_blank(),

    axis.text.x.bottom = element_text(size = 9, color = "black"),
    axis.title.x.bottom = element_text(size = 9, color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x.bottom = element_line(color = "black")
  ) +

  scale_y_continuous(position = "left") +
  scale_x_discrete(position = "top")


ggsave(file.path(existing_results_path, "Figure2C_Barplot.pdf"), panel1, width = 6, height = 4)


survival_data$DiffCat <- cut(abs(survival_data$Difference),
                             breaks = c(-Inf, 0.03, 0.05, 0.08, Inf),
                             labels = c("< 0.03", "0.03-0.05", "0.05-0.08", "> 0.08"))

survival_data$Label <- ifelse(abs(survival_data$Difference) > 0.08, survival_data$State, "")


scatterplot_data <- data.frame(
  State = survival_data$State,
  Alive_Score = survival_data$Alive,
  Dead_Score = survival_data$Dead,
  Difference = survival_data$Difference,
  Difference_Category = survival_data$DiffCat,
  Label = survival_data$Label
)
write.csv(scatterplot_data, file.path(existing_results_path, "Figure2D_Scatterplot_Data.csv"), row.names = FALSE)
log_message("Scatter plot data table saved\n", print_to_console = TRUE)

panel2 <- ggplot(survival_data, aes(x = Alive, y = Dead, color = Difference, size = abs(Difference))) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = Label), size = 2.5, vjust = -0.8, hjust = 0.5, check_overlap = TRUE, color = "black") +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 4)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgray") +
  theme_light() +
  labs(
    title = "Functional State Activity: Alive vs Dead",
    x = "Alive Group Average Activity",
    y = "Dead Group Average Activity",
    color = "Difference",
    size = "Absolute Difference"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 9, color = "black"),
    legend.text = element_text(size = 8, color = "black")
  )


ggsave(file.path(existing_results_path, "Figure2D_Scatterplot.pdf"), panel2, width = 8, height = 6)


log_message("Creating volcano plot...
", print_to_console = TRUE)


volcano_data <- survival_data
volcano_data$logFC <- volcano_data$Difference
volcano_data$negLogP <- abs(volcano_data$Difference) * 20


volcano_data$Regulation <- ifelse(volcano_data$Difference > 0.05, "Upregulated",
                                 ifelse(volcano_data$Difference < -0.05, "Downregulated", "No Change"))


volcano_data$Label <- ifelse(abs(volcano_data$Difference) > 0.08, volcano_data$State, "")


volcano_data_output <- data.frame(
  State = volcano_data$State,
  logFC = volcano_data$logFC,
  negLogP = volcano_data$negLogP,
  Regulation = volcano_data$Regulation,
  Label = volcano_data$Label
)
write.csv(volcano_data_output, file.path(existing_results_path, "Figure2E_Volcano_Data.csv"), row.names = FALSE)
log_message("Volcano plot data table saved\n", print_to_console = TRUE)


if (!requireNamespace("ggrepel", quietly = TRUE)) {

  panel3 <- ggplot(volcano_data, aes(x = logFC, y = negLogP, color = Regulation)) +
    geom_point(aes(size = abs(Difference)), alpha = 0.7) +
    geom_text(aes(label = Label), size = 4, vjust = -0.5, hjust = 0.5) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No Change" = "gray")) +
    theme_light() +
    labs(
      title = "Differences in Functional States Between Dead and Alive Groups",
      x = "Fold Change (Dead - Alive)",
      y = "Significance (-log10 p-value simulation)",
      color = "Regulation"
    ) +
    geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "darkgray") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
    theme(
      plot.title = element_text(size = 12, face = "bold", color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      axis.text = element_text(size = 8, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 9, color = "black"),
      legend.text = element_text(size = 8, color = "black")
    )
  } else {

    library(ggrepel)
    panel3 <- ggplot(volcano_data, aes(x = logFC, y = negLogP, color = Regulation, label = Label)) +
      geom_point(aes(size = abs(Difference)), alpha = 0.7) +
      geom_text_repel(size = 4, max.overlaps = 15, box.padding = 0.5) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No Change" = "gray")) +
      theme_light() +
      labs(
        title = "Differences in Functional States Between Dead and Alive Groups",
        x = "Fold Change (Dead - Alive)",
        y = "Significance (-log10 p-value simulation)",
        color = "Regulation"
      ) +
      geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "darkgray") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
      theme(
        plot.title = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 8, color = "black")
      )
  }


ggsave(file.path(existing_results_path, "Figure2E_Volcano_Plot.pdf"), panel3, width = 8, height = 6)
log_message("Volcano plot saved\n", print_to_console = TRUE)


if (!is.null(cox_results)) {
  top_cox_states <- intersect(top_states$State, cox_results$State)
  cox_subset <- cox_results[cox_results$State %in% top_cox_states, ]


  cox_subset$State <- factor(cox_subset$State,
                            levels = rev(cox_subset$State[order(cox_subset$HR)]))
  cox_subset$Significance <- ifelse(cox_subset$P_Value < 0.05, "Significant", "Non-significant")


  cox_forest_data <- data.frame(
    State = levels(cox_subset$State),
    HR = cox_subset$HR[match(levels(cox_subset$State), cox_subset$State)],
    CI_Lower = cox_subset$CI_Lower[match(levels(cox_subset$State), cox_subset$State)],
    CI_Upper = cox_subset$CI_Upper[match(levels(cox_subset$State), cox_subset$State)],
    P_Value = cox_subset$P_Value[match(levels(cox_subset$State), cox_subset$State)],
    Significance = cox_subset$Significance[match(levels(cox_subset$State), cox_subset$State)]
  )
  write.csv(cox_forest_data, file.path(existing_results_path, "Figure2F_Cox_Forest_Data.csv"), row.names = FALSE)
  log_message("Cox regression forest plot data table saved, order consistent with figure\n", print_to_console = TRUE)

  panel4 <- ggplot(cox_subset, aes(x = HR, y = State, xmin = CI_Lower, xmax = CI_Upper, color = Significance)) +
    geom_point(size = 3) +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray50")) +
    theme_light() +
    labs(
      title = "Cox Regression Analysis (Hazard Ratio)",
      x = "Hazard Ratio (HR)",
      y = "Functional State",
      color = "Statistical Significance"
    ) +
    theme(
      plot.title = element_text(size = 9, face = "bold", color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "top",
      legend.title = element_text(size = 9, color = "black"),
      legend.text = element_text(size = 9, color = "black"),

      axis.text.y.right = element_text(size = 9, color = "black"),
      axis.title.y.right = element_text(size = 9, color = "black", angle = 90),

      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),

      axis.line.y.right = element_line(color = "black"),
      axis.ticks.y.right = element_line(color = "black")
    ) +

    scale_y_discrete(position = "right")


  ggsave(file.path(existing_results_path, "Figure2F_Cox_Forest_Plot.pdf"), panel4, width = 5.7, height = 4)
}


top_states_data <- data.frame(
  State = as.character(top_states$State),
  Difference = top_states$Difference[match(levels(top_states$State), top_states$State)],
  Direction = top_states$Direction[match(levels(top_states$State), top_states$State)],
  Alive_Score = top_states$Alive[match(levels(top_states$State), top_states$State)],
  Dead_Score = top_states$Dead[match(levels(top_states$State), top_states$State)],
  AbsDiff = top_states$AbsDiff[match(levels(top_states$State), top_states$State)]
)
write.csv(top_states_data, file.path(existing_results_path, "Figure2C_Top10_States.csv"), row.names = FALSE)
log_message("Top 10 functional state data table saved, order consistent with figure\n", print_to_console = TRUE)


log_message("\nAll figures and data files saved separately!
", print_to_console = TRUE)
log_message("Saved files:
", print_to_console = TRUE)
log_message("   - Figure2C_Top10_States.csv\n", print_to_console = TRUE)
log_message("   - Figure2D_Scatterplot_Data.csv\n", print_to_console = TRUE)
log_message("   - Figure2F_Cox_Forest_Data.csv\n", print_to_console = TRUE)

log_message("Analysis complete!
", print_to_console = TRUE)
