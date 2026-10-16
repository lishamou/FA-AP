

log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "Figure4_key_comparison.log")
log_con <- file(log_file, "w")


log_message <- function(message, print_to_console = FALSE) {
  cat(message, file = log_con, append = TRUE)
  if (print_to_console) {
    cat(message)
  }
}


library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)


output_dir <- "./results"


risk_group_comparison <- read.csv(file.path(output_dir, "Figure4_Risk_Group_Comparison.csv"))
risk_group_long_format <- read.csv(file.path(output_dir, "Figure4_Risk_Group_Long_Format.csv"))


risk_group_comparison <- risk_group_comparison %>%
  arrange(P_Value)


top_states <- head(risk_group_comparison$AgentState, 6)
cat("The 6 most significant agent states:\n", paste(top_states, collapse = "\n"), "\n")


high_risk_enriched <- c(
  "CSC_Activity_PRG",
  "Autophagy_Endothelial_HighRisk_AutRG",
  "CytokineSignaling_HighRisk_CSIRG",
  "KEGG_HOMOLOGOUS_RECOMBINATION",
  "KEGG_CELL_CYCLE",
  "KEGG_DNA_REPLICATION"
)


low_risk_enriched <- c(
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "LipidMetabolism_Risk_LMRG",
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS",
  "Mast_cells",
  "B_cells"
)


key_states <- c(high_risk_enriched, low_risk_enriched)


key_states <- unique(key_states)


filtered_data <- risk_group_long_format %>%
  filter(AgentState %in% key_states)


create_comparison_plot <- function(data, state_name) {

  state_data <- data[data$AgentState == state_name, ]


  p_val <- tryCatch({
    wilcox.test(Score ~ RiskGroup, data = state_data)$p.value
  }, error = function(e) {
    cat("Warning: Unable to calculate p-value for", state_name, ":", e$message, "\n")
    return(NA)
  })


  if (is.na(p_val)) {
    p_text <- "p = NA"
  } else if (p_val < 0.001) {
    p_text <- "p < 0.001"
  } else if (p_val < 0.01) {
    p_text <- paste0("p = ", sprintf("%.3f", p_val))
  } else {
    p_text <- paste0("p = ", sprintf("%.3f", p_val))
  }


  y_max <- max(state_data$Score, na.rm = TRUE)
  y_min <- min(state_data$Score, na.rm = TRUE)
  y_range <- y_max - y_min
  y_pos <- y_max + 0.1 * y_range


  display_name <- state_name


  high_risk_colors <- c("Low" = "#2C7BB6", "High" = "#D7191C")
  low_risk_colors <- c("Low" = "#9ab55e", "High" = "#bd7cb3")


  enriched_in <- ifelse(state_name %in% high_risk_enriched, "High Risk", "Low Risk")


  if (state_name %in% high_risk_enriched) {
    risk_colors <- high_risk_colors
  } else {
    risk_colors <- low_risk_colors
  }


  p <- ggplot(state_data, aes(x = RiskGroup, y = Score, fill = RiskGroup)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = risk_colors) +

    scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
    labs(
      title = display_name,
      x = NULL,
      y = "Activity Score (GSVA)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",

      text = element_text(color = "black"),

      plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 12, color = "black"),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(colour = "black", fill = NA)
    ) +
    scale_x_discrete(labels = c("Low" = "Low Risk", "High" = "High Risk")) +
    annotate("text", x = 1.5, y = y_pos, label = p_text, size = 5, color = "black")

  return(p)
}


plot_list <- list()
for (state in key_states) {
  cat("Drawing plot for", state, "...\n")
  plot_list[[state]] <- create_comparison_plot(filtered_data, state)
}


high_risk_plots <- plot_list[high_risk_enriched]
low_risk_plots <- plot_list[low_risk_enriched]


cat("Saving plots for high risk enriched states...\n")


high_risk_data <- filtered_data %>%
  filter(AgentState %in% high_risk_enriched) %>%
  select(AgentState, RiskGroup, Score, Cohort) %>%
  arrange(AgentState, RiskGroup)
write.csv(high_risk_data, file.path(output_dir, "Figure4A_HighRisk_Enriched_Data.csv"), row.names = FALSE)
cat("Saved high risk enriched data to CSV\n")


high_risk_pdf <- file.path(output_dir, "Figure4A_HighRisk_Enriched_Combined.pdf")
pdf(high_risk_pdf, width = 12, height = 8)
grid.arrange(grobs = high_risk_plots, ncol = 3)
dev.off()


cat("Saving plots for low risk enriched states...\n")


low_risk_data <- filtered_data %>%
  filter(AgentState %in% low_risk_enriched) %>%
  select(AgentState, RiskGroup, Score, Cohort) %>%
  arrange(AgentState, RiskGroup)
write.csv(low_risk_data, file.path(output_dir, "Figure4B_LowRisk_Enriched_Data.csv"), row.names = FALSE)
cat("Saved low risk enriched data to CSV\n")


low_risk_pdf <- file.path(output_dir, "Figure4B_LowRisk_Enriched_Combined.pdf")
pdf(low_risk_pdf, width = 12, height = 8)
grid.arrange(grobs = low_risk_plots, ncol = 3)
dev.off()


cat("Combined plots have been saved to", output_dir, "directory\n")
