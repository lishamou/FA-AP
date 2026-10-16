
library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(stringr)


source("code/theme_sci_blue.R")


correlation_matrix <- read.csv("results/survival_prediction_optimized_py/2_2.cleaned_feature_correlation.csv",
                             row.names = 1, check.names = FALSE)


simplify_feature_names <- function(names) {
  short_names <- names %>%
    str_replace("CytokineSignaling_HighRisk_CSIRG", "CSIRG") %>%
    str_replace("Metabolism_General_Risk_MRG", "MRG") %>%
    str_replace("Autophagy_Endothelial_HighRisk_AutRG", "AutRG") %>%
    str_replace("CSC_Activity_PRG", "PRG") %>%
    str_replace("KEGG_", "K_") %>%
    str_replace("GLYCINE_SERINE_AND_THREONINE_METABOLISM", "GLY_SER_THR") %>%
    str_replace("_METABOLISM", "_MET") %>%
    str_replace("_SIGNALING_PATHWAY", "_SIG") %>%
    str_replace("_REPLICATION", "_REP") %>%
    str_replace("REGULATION_OF_", "REG_") %>%
    str_replace("Type_I_IFN_Reponse", "IFN-I") %>%
    str_replace("Type_II_IFN_Reponse", "IFN-II") %>%
    str_replace("COMPLEMENT_AND_COAGULATION_CASCADES", "COMP_COAG")
  return(short_names)
}


simplified_names <- simplify_feature_names(colnames(correlation_matrix))


name_mapping <- data.frame(
  Original = colnames(correlation_matrix),
  Simplified = simplified_names,
  stringsAsFactors = FALSE
)


colnames(correlation_matrix) <- rownames(correlation_matrix) <- simplified_names


output_dir <- "results/survival_prediction_optimized_py"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(name_mapping, file.path(output_dir, "2_2.cleaned_feature_name_mapping.csv"), row.names = FALSE)
cat("Feature name mapping table saved to: ", file.path(output_dir, "2_2.cleaned_feature_name_mapping.csv"), "\n")


melted_corr_matrix <- melt(as.matrix(correlation_matrix))
colnames(melted_corr_matrix) <- c("Feature1", "Feature2", "Correlation")


p_full <- ggplot(data = melted_corr_matrix, aes(x = Feature1, y = Feature2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +

  scale_fill_gradient2(
    low = "#053061",
    mid = "#f7f7f7",
    high = "#67001f",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.2),
    name = "Correlation Coefficient",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barwidth = 1,
      barheight = 10
    )
  ) +

  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Feature Correlation Matrix (Full Matrix)",
    x = "",
    y = ""
  ) +

  scale_x_discrete(limits = colnames(correlation_matrix)) +
  scale_y_discrete(limits = rev(colnames(correlation_matrix)))


lower_tri <- melted_corr_matrix %>%
  filter(as.numeric(factor(Feature1, levels = colnames(correlation_matrix))) <=
         as.numeric(factor(Feature2, levels = colnames(correlation_matrix))))


common_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )


p1 <- ggplot(data = lower_tri, aes(x = Feature1, y = Feature2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient2(
    low = "#053061",
    mid = "#f7f7f7",
    high = "#67001f",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation Coefficient\n(Red-Blue)"
  ) +
  common_theme +
  labs(
    title = "Feature Correlation Matrix (Red-Blue)",
    x = "",
    y = ""
  ) +
  scale_x_discrete(limits = colnames(correlation_matrix)) +
  scale_y_discrete(limits = rev(colnames(correlation_matrix)))


p2 <- ggplot(data = lower_tri, aes(x = Feature1, y = Feature2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient2(
    low = "#5e4fa2",
    mid = "#f7f7f7",
    high = "#d73027",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation Coefficient\n(Purple-Orange)"
  ) +
  common_theme +
  labs(
    title = "Feature Correlation Matrix (Purple-Orange)",
    x = "",
    y = ""
  ) +
  scale_x_discrete(limits = colnames(correlation_matrix)) +
  scale_y_discrete(limits = rev(colnames(correlation_matrix)))


p3 <- ggplot(data = lower_tri, aes(x = Feature1, y = Feature2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_distiller(
    palette = "BrBG",
    limits = c(-1, 1),
    name = "Correlation Coefficient\n(BrBG)"
  ) +
  common_theme +
  labs(
    title = "Feature Correlation Matrix (BrBG)",
    x = "",
    y = ""
  ) +
  scale_x_discrete(limits = colnames(correlation_matrix)) +
  scale_y_discrete(limits = rev(colnames(correlation_matrix)))


p4 <- ggplot(data = lower_tri, aes(x = Feature1, y = Feature2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradientn(
    colors = sci_blue_palette(100),
    limits = c(-1, 1),
    name = "Correlation Coefficient\n(SCI Blue)"
  ) +
  sci_theme_blue() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Feature Correlation Matrix (SCI Blue)",
    x = "",
    y = ""
  ) +
  scale_x_discrete(limits = colnames(correlation_matrix)) +
  scale_y_discrete(limits = rev(colnames(correlation_matrix)))


ggsave(file.path(output_dir, "Feature_Correlation_Matrix_LowerTri_RedBlue.pdf"), p1, width = 14, height = 12)

ggsave(file.path(output_dir, "Feature_Correlation_Matrix_LowerTri_PurpleOrange.pdf"), p2, width = 14, height = 12)

ggsave(file.path(output_dir, "Feature_Correlation_Matrix_LowerTri_BrBG.pdf"), p3, width = 14, height = 12)


ggsave(file.path(output_dir, "Figure6A_Correlation_Matrix_LowerTri_SCI_Blue.pdf"), p4, width = 14, height = 12)


cat("\nFigure saved to:
")
cat("2. PDF format: results/survival_prediction_optimized_py/Feature_Correlation_Matrix_R.pdf\n")
