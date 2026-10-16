
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)


source("code/theme_sci_blue.R")


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

output_dir <- "results/survival_prediction_optimized_py"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


feature_data <- read.csv(file.path(output_dir, "2_2.Cleaned_Feature_Importance.csv"))


feature_data$Feature_Short <- simplify_feature_names(feature_data$Feature)


cat('Original feature name: ', feature_data$Feature, '\n')
cat('Simplified feature name: ', feature_data$Feature_Short, '\n')
cat('Length of simplified feature name: ', length(feature_data$Feature_Short), '\n')
cat('Number of data rows: ', nrow(feature_data), '\n')


feature_data$Feature_Short <- as.character(feature_data$Feature_Short)


feature_data <- feature_data %>% arrange(desc(Importance))


complete_data_file <- file.path(output_dir, "2_2.Cleaned_Feature_Importance_with_ShortNames.csv")
write.csv(feature_data, complete_data_file, row.names = FALSE)
cat("Complete feature data (with simplified names) saved to: ", complete_data_file, "\n")


name_mapping <- data.frame(Original = feature_data$Feature, Simplified = feature_data$Feature_Short, stringsAsFactors = FALSE)
write.csv(name_mapping, file.path(output_dir, "2_2.cleaned_feature_name_mapping.csv"), row.names = FALSE)
cat("Feature name mapping table saved to: ", file.path(output_dir, "2_2.cleaned_feature_name_mapping.csv"), "\n")


p <- ggplot(feature_data, aes(x = reorder(Feature_Short, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = sci_blue_palette(1)[1]) +
  coord_flip() +
  sci_theme_blue() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "Importance Distribution of Final Selected Features",
    x = "",
    y = "Feature Importance"
  )

print(p)


ggsave(
  file.path(output_dir, "2_2.Cleaned_Feature_Importance.pdf"),
  plot = p,
  width = 4,
  height = 8
)

cat("\nFeature importance ranking:
")
print(feature_data %>% select(Feature_Short, Importance) %>% mutate(Importance = round(Importance, 4)))

cat("\nFigure and data saved to:
")
cat("1. Complete data (with simplified names): ", complete_data_file, "\n")
cat("2. Feature name mapping table: ", file.path(output_dir, "2_2.cleaned_feature_name_mapping.csv"), "\n")
cat("3. PDF-format figure: ", file.path(output_dir, "2_2.Cleaned_Feature_Importance.pdf"), "\n")


feature_data <- read.csv(file.path(output_dir, "4.cv_selected_features.csv"))
feature_data$Feature_Short <- simplify_feature_names(feature_data$Feature)


feature_data <- feature_data %>% arrange(desc(Importance))


complete_data_file <- file.path(output_dir, "4.cv_selected_features_with_ShortNames.csv")
write.csv(feature_data, complete_data_file, row.names = FALSE)
cat("\nCV-selected complete data (with simplified names) saved to: ", complete_data_file, "\n")


name_mapping <- data.frame(Original = feature_data$Feature, Simplified = feature_data$Feature_Short, stringsAsFactors = FALSE)
write.csv(name_mapping, file.path(output_dir, "4.cv_selected_features_name_mapping.csv"), row.names = FALSE)
cat("CV-selected name mapping table saved to: ", file.path(output_dir, "4.cv_selected_features_name_mapping.csv"), "\n")

p <- ggplot(feature_data, aes(x = reorder(Feature_Short, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = sci_blue_palette(1)[1]) +
  coord_flip() +
  sci_theme_blue() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "Importance Distribution of CV Selected Features",
    x = "",
    y = "Feature Importance"
  )

print(p)

ggsave(
  file.path(output_dir, "Figure6C_cv_selected_features.pdf"),
  plot = p,
  width = 4,
  height = 5
)

cat("\nFeature importance ranking:
")
print(feature_data %>% select(Feature_Short, Importance) %>% mutate(Importance = round(Importance, 4)))

cat("\nCV-selected figure and data saved to:
")
cat("1. Complete data (with simplified names): ", complete_data_file, "\n")
cat("2. Feature name mapping table: ", file.path(output_dir, "4.cv_selected_features_name_mapping.csv"), "\n")
cat("3. PDF-format figure: ", file.path(output_dir, "Figure6C_cv_selected_features.pdf"), "\n")


feature_data <- read.csv(file.path(output_dir, "5.final_feature_list.csv"))
feature_data$Feature_Short <- simplify_feature_names(feature_data$Feature)

feature_data <- feature_data %>% arrange(desc(Importance))


complete_data_file <- file.path(output_dir, "5.final_feature_list_with_ShortNames.csv")
write.csv(feature_data, complete_data_file, row.names = FALSE)
cat("\nFinal selected complete data (with simplified names) saved to: ", complete_data_file, "\n")


name_mapping <- data.frame(Original = feature_data$Feature, Simplified = feature_data$Feature_Short, stringsAsFactors = FALSE)
write.csv(name_mapping, file.path(output_dir, "5.final_feature_list_name_mapping.csv"), row.names = FALSE)
cat("Final selected name mapping table saved to: ", file.path(output_dir, "5.final_feature_list_name_mapping.csv"), "\n")

p <- ggplot(feature_data, aes(x = reorder(Feature_Short, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = sci_blue_palette(1)[1]) +
  coord_flip() +
  sci_theme_blue() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "Importance of Final Selected Features",
    x = "",
    y = "Feature Importance"
  )

print(p)

ggsave(
  file.path(output_dir, "Figure6D_final_feature_list.pdf"),
  plot = p,
  width = 5,
  height = 8
)

cat("\nFeature importance ranking:
")
print(feature_data %>% select(Feature_Short, Importance) %>% mutate(Importance = round(Importance, 4)))

cat("\nFinal selected figure and data saved to:
")
cat("1. Complete data (with simplified names): ", complete_data_file, "\n")
cat("2. Feature name mapping table: ", file.path(output_dir, "5.final_feature_list_name_mapping.csv"), "\n")
cat("3. PDF-format figure: ", file.path(output_dir, "Figure6D_final_feature_list.pdf"), "\n")
