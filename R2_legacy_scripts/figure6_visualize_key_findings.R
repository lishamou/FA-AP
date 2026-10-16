

Sys.setlocale("LC_ALL", "en_US.UTF-8")


log_dir <- "log"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
log_file <- file.path(log_dir, "Figure6_visualize_key_findings.log")
log_con <- file(log_file, "w")


log_message <- function(message, print_to_console = TRUE) {
  cat(message, file = log_con, append = TRUE)
  if (print_to_console) {
    cat(message)
  }
}


on.exit({
  close(log_con)
  log_message(paste0("\nLog saved to: ", normalizePath(log_file), "\n"), print_to_console = TRUE)
})


load_packages <- function() {
  log_message("Loading required R packages...
", print_to_console = TRUE)

  required_packages <- c("ggplot2", "dplyr", "tidyr", "ggpubr", "ComplexHeatmap",
                         "circlize", "RColorBrewer", "patchwork", "viridis", "ggrepel")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      log_message(paste0("Installing package: ", pkg, "\n"), print_to_console = TRUE)
      install.packages(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }


  dplyr_version <- as.character(packageVersion("dplyr"))
  log_message(paste0("dplyr version used: ", dplyr_version, "\n"), print_to_console = TRUE)


  tidyr_version <- as.character(packageVersion("tidyr"))
  log_message(paste0("tidyr version used: ", tidyr_version, "\n"), print_to_console = TRUE)

  log_message("All required R packages loaded
", print_to_console = TRUE)
}


theme_publication <- function(base_size = 18, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      axis.text = element_text(colour = "black", size = base_size * 0.8),
      axis.title = element_text(size = base_size),
      legend.title = element_text(size = base_size * 0.8),
      legend.text = element_text(size = base_size * 0.7),
      legend.position = "right",
      legend.key.size = unit(0.8, "lines"),
      plot.title = element_text(size = base_size * 1.2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size * 0.8, hjust = 0.5),
      plot.caption = element_text(size = base_size * 0.6, hjust = 1)
    )
}


read_data <- function() {
  log_message("Reading CSV result file...
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
        file.exists(file.path(path, "Figure6_TP53_Mutation_Comparison.csv")) &&
        file.exists(file.path(path, "Figure6_CTNNB1_Mutation_Comparison.csv"))) {
      existing_results_path <- path
      log_message(paste0("Found data file directory: ", existing_results_path, "\n"), print_to_console = TRUE)
      break
    }
  }


  if (is.null(existing_results_path)) {
    stop("Error: previously generated data file not found. Please run Figure6.R first.")
  }


  if (!dir.exists(existing_results_path)) {
    dir.create(existing_results_path, recursive = TRUE)
    log_message(paste0("Creating output directory: ", existing_results_path, "\n"), print_to_console = TRUE)
  }


  tp53_data <- read.csv(file.path(existing_results_path, "Figure6_TP53_Mutation_Comparison.csv"))
  ctnnb1_data <- read.csv(file.path(existing_results_path, "Figure6_CTNNB1_Mutation_Comparison.csv"))

  log_message(paste0("Successfully read TP53 data, ", nrow(tp53_data), " rows
"), print_to_console = TRUE)
  log_message(paste0("Successfully read CTNNB1 data, ", nrow(ctnnb1_data), " rows
"), print_to_console = TRUE)

  return(list(tp53_data = tp53_data, ctnnb1_data = ctnnb1_data, output_dir = existing_results_path))
}


filter_significant_results <- function(data, p_value_threshold = 0.05, log2fc_threshold = 0.005) {
  log_message(paste0("Filtering P<", p_value_threshold, " and |Log2FC|>", log2fc_threshold, " significant results...\n"),
              print_to_console = TRUE)


  data$Change <- data$Mean_Score_Mutant - data$Mean_Score_WT
  data$Log2FC <- data$Change


  significant_data <- tryCatch({

    data %>%
      filter(!is.na(P_Value) & P_Value < p_value_threshold & abs(Change) > log2fc_threshold) %>%
      arrange(P_Value)
  }, error = function(e) {
    log_message(paste0("dplyr filter/arrange failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


    filtered_data <- data[!is.na(data$P_Value) &
                          data$P_Value < p_value_threshold &
                          abs(data$Change) > log2fc_threshold, ]
    filtered_data <- filtered_data[order(filtered_data$P_Value), ]
    return(filtered_data)
  })

  log_message(paste0("Filtered ", nrow(significant_data), " significant differential states\n"), print_to_console = TRUE)

  return(significant_data)
}


plot_volcano <- function(data, gene_name, output_dir) {
  log_message(paste0("Plotting ", gene_name, " mutation-related volcano plot...\n"), print_to_console = TRUE)


  plot_data <- tryCatch({

    data %>%
      mutate(
        Change = Mean_Score_Mutant - Mean_Score_WT,
        neg_log10_p = -log10(P_Value),
        Significance = case_when(
          is.na(P_Value) ~ "NS",
          P_Value < 0.05 & Change > 0.01 ~ "Up",
          P_Value < 0.05 & Change < -0.01 ~ "Down",
          TRUE ~ "NS"
        )
      )
  }, error = function(e) {
    log_message(paste0("dplyr mutate/case_when failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


    result <- data
    result$Change <- result$Mean_Score_Mutant - result$Mean_Score_WT
    result$neg_log10_p <- -log10(result$P_Value)


    result$Significance <- "NS"
    result$Significance[!is.na(result$P_Value) & result$P_Value < 0.05 & result$Change > 0.01] <- "Up"
    result$Significance[!is.na(result$P_Value) & result$P_Value < 0.05 & result$Change < -0.01] <- "Down"

    return(result)
  })


  top_up <- plot_data %>%
    filter(Significance == "Up") %>%
    arrange(P_Value) %>%
    head(5)

  top_down <- plot_data %>%
    filter(Significance == "Down") %>%
    arrange(P_Value) %>%
    head(5)


  top_results <- bind_rows(top_up, top_down)
  plot_data$Label <- ifelse(plot_data$AgentState %in% top_results$AgentState,
                           as.character(plot_data$AgentState),
                           "")


  volcano_colors <- c(
    "Up" = "#D55E00",
    "Down" = "#0072B2",
    "NS" = "grey80"
  )


  p <- ggplot(plot_data, aes(x = Change, y = neg_log10_p, color = Significance)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = volcano_colors) +
    ggrepel::geom_text_repel(
      data = plot_data %>% filter(Label != ""),
      aes(label = Label),
      size = 5,
      max.overlaps = 10,
      box.padding = 0.5,
      segment.alpha = 0.5,
      min.segment.length = 0,
      force = 10
    ) +
    labs(
      title = paste0(gene_name, " Mutation-Related Agent State Changes"),
      subtitle = "Top 5 up- and down-regulated states labeled",
      x = "Activity Score Change (Mutant - Wild Type)",
      y = "-Log10(P-value)",
      color = "Change Direction"
    ) +
    theme_publication() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )


  volcano_csv <- file.path(output_dir, paste0("Figure6AB_", gene_name, "_Volcano_Data.csv"))
  write.csv(plot_data, volcano_csv, row.names = FALSE)
  log_message(paste0("Volcano plot data saved: ", volcano_csv, "\n"), print_to_console = TRUE)


  top_results_csv <- file.path(output_dir, paste0("Figure6_", gene_name, "_Volcano_Top5_Results.csv"))
  write.csv(top_results, top_results_csv, row.names = FALSE)
  log_message(paste0("Top 5 up/down results saved: ", top_results_csv, "\n"), print_to_console = TRUE)

  return(p)
}


plot_heatmap <- function(tp53_significant, ctnnb1_significant, output_dir, top_n = 100,
                         row_fontsize = 10, col_fontsize = 12,
                         legend_title_fontsize = 10, legend_labels_fontsize = 8,
                         star_fontsize = 10) {
  log_message(paste0("Plotting heatmap of top ", top_n, " TP53/CTNNB1 mutation-related states..."), print_to_console = TRUE)


  log_message(paste0("TP53 data row count: ", nrow(tp53_significant), "\n"), print_to_console = TRUE)
  log_message(paste0("CTNNB1 data row count: ", nrow(ctnnb1_significant), "\n"), print_to_console = TRUE)


  required_columns <- c("AgentState", "Mean_Score_Mutant", "Mean_Score_WT", "P_Value")
  for (col in required_columns) {
    if (!(col %in% colnames(tp53_significant)) || !(col %in% colnames(ctnnb1_significant))) {
      error_msg <- paste0("Error: required column '", col, "'")
      log_message(error_msg, print_to_console = TRUE)
      stop(error_msg)
    }
  }


  if (!("Change" %in% colnames(tp53_significant))) {
    tp53_significant$Change <- tp53_significant$Mean_Score_Mutant - tp53_significant$Mean_Score_WT
    log_message("Adding Change column to TP53 data\n", print_to_console = TRUE)
  }
  if (!("Change" %in% colnames(ctnnb1_significant))) {
    ctnnb1_significant$Change <- ctnnb1_significant$Mean_Score_Mutant - ctnnb1_significant$Mean_Score_WT
    log_message("Adding Change column to CTNNB1 data\n", print_to_console = TRUE)
  }


  top_tp53 <- tryCatch({

    actual_top_n <- min(top_n, nrow(tp53_significant))
    log_message(paste0("Actual TP53 top_n used: ", actual_top_n, "\n"), print_to_console = TRUE)


    tp53_significant %>%
      arrange(P_Value) %>%
      head(actual_top_n)
  }, error = function(e) {
    log_message(paste0("dplyr arrange/head failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


    ordered_indices <- order(tp53_significant$P_Value)
    actual_top_n <- min(top_n, length(ordered_indices))
    tp53_significant[ordered_indices[1:actual_top_n], ]
  })

  top_ctnnb1 <- tryCatch({

    actual_top_n <- min(top_n, nrow(ctnnb1_significant))
    log_message(paste0("Actual CTNNB1 top_n used: ", actual_top_n, "\n"), print_to_console = TRUE)


    ctnnb1_significant %>%
      arrange(P_Value) %>%
      head(actual_top_n)
  }, error = function(e) {
    log_message(paste0("dplyr arrange/head failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


    ordered_indices <- order(ctnnb1_significant$P_Value)
    actual_top_n <- min(top_n, length(ordered_indices))
    ctnnb1_significant[ordered_indices[1:actual_top_n], ]
  })


  all_states <- unique(c(top_tp53$AgentState, top_ctnnb1$AgentState))
  log_message(paste0("Number of merged states: ", length(all_states), "\n"), print_to_console = TRUE)
  log_message("State list:
", print_to_console = TRUE)
  for (state in all_states) {
    log_message(paste0(" - ", state, "\n"), print_to_console = TRUE)
  }


  heatmap_data <- data.frame(
    AgentState = all_states,
    stringsAsFactors = FALSE
  )


  for (state in all_states) {

    tp53_row <- tryCatch({

      filtered_data <- tp53_significant %>% filter(AgentState == state)
      log_message(paste0("Processing TP53 state '", state, "', Found ", nrow(filtered_data), " rows
"), print_to_console = TRUE)
      filtered_data
    }, error = function(e) {

      filtered_data <- tp53_significant[tp53_significant$AgentState == state, ]
      log_message(paste0("Processing TP53 state '", state, "', Found ", nrow(filtered_data), " rows
"), print_to_console = TRUE)
      filtered_data
    })

    if (nrow(tp53_row) > 0) {
      heatmap_data[heatmap_data$AgentState == state, "TP53_Change"] <- tp53_row$Change[1]
      heatmap_data[heatmap_data$AgentState == state, "TP53_P_Value"] <- tp53_row$P_Value[1]
    } else {
      heatmap_data[heatmap_data$AgentState == state, "TP53_Change"] <- 0
      heatmap_data[heatmap_data$AgentState == state, "TP53_P_Value"] <- 1
      log_message(paste0("Warning: TP53 state '", state, "'' not found, using default value
"), print_to_console = TRUE)
    }


    ctnnb1_row <- tryCatch({

      filtered_data <- ctnnb1_significant %>% filter(AgentState == state)
      log_message(paste0("Processing CTNNB1 state '", state, "', Found ", nrow(filtered_data), " rows
"), print_to_console = TRUE)
      filtered_data
    }, error = function(e) {

      filtered_data <- ctnnb1_significant[ctnnb1_significant$AgentState == state, ]
      log_message(paste0("Processing CTNNB1 state '", state, "', Found ", nrow(filtered_data), " rows
"), print_to_console = TRUE)
      filtered_data
    })

    if (nrow(ctnnb1_row) > 0) {
      heatmap_data[heatmap_data$AgentState == state, "CTNNB1_Change"] <- ctnnb1_row$Change[1]
      heatmap_data[heatmap_data$AgentState == state, "CTNNB1_P_Value"] <- ctnnb1_row$P_Value[1]
    } else {
      heatmap_data[heatmap_data$AgentState == state, "CTNNB1_Change"] <- 0
      heatmap_data[heatmap_data$AgentState == state, "CTNNB1_P_Value"] <- 1
      log_message(paste0("Warning: CTNNB1 state '", state, "'' not found, using default value
"), print_to_console = TRUE)
    }
  }


  log_message("\nHeatmap data summary:
", print_to_console = TRUE)
  log_message(paste0("Rows: ", nrow(heatmap_data), "\n"), print_to_console = TRUE)
  log_message(paste0("Columns: ", ncol(heatmap_data), "\n"), print_to_console = TRUE)
  log_message(paste0("Column names: ", toString(colnames(heatmap_data)), "\n"), print_to_console = TRUE)


  required_cols <- c("AgentState", "TP53_Change", "CTNNB1_Change", "TP53_P_Value", "CTNNB1_P_Value")
  missing_cols <- setdiff(required_cols, colnames(heatmap_data))
  if (length(missing_cols) > 0) {
    error_msg <- paste0("Error: required column missing in heatmap data: ", paste(missing_cols, collapse = ", "))
    log_message(paste0(error_msg, "\n"), print_to_console = TRUE)
    stop(error_msg)
  }


  heatmap_matrix <- as.matrix(heatmap_data[, c("TP53_Change", "CTNNB1_Change")])
  rownames(heatmap_matrix) <- heatmap_data$AgentState


  p_value_matrix <- as.matrix(heatmap_data[, c("TP53_P_Value", "CTNNB1_P_Value")])
  rownames(p_value_matrix) <- heatmap_data$AgentState


  log_message("\nMatrix dimensions:
", print_to_console = TRUE)
  log_message(paste0("Heatmap matrix: ", paste(dim(heatmap_matrix), collapse = " x "), "\n"), print_to_console = TRUE)
  log_message(paste0("P-value matrix: ", paste(dim(p_value_matrix), collapse = " x "), "\n"), print_to_console = TRUE)


  if (any(is.na(heatmap_matrix))) {
    log_message("Warning: heatmap matrix contains NA values
", print_to_console = TRUE)
  }
  if (any(is.na(p_value_matrix))) {
    log_message("Warning: P-value matrix contains NA values
", print_to_console = TRUE)
  }


  avg_change <- rowMeans(heatmap_matrix)

  sort_order <- order(
    sign(avg_change),
    -abs(avg_change)
  )


  heatmap_matrix <- heatmap_matrix[sort_order, ]
  p_value_matrix <- p_value_matrix[sort_order, ]


  max_abs_value <- max(abs(heatmap_matrix))
  col_fun <- colorRamp2(
    c(-max_abs_value, -max_abs_value/2, -max_abs_value/4, 0, max_abs_value/4, max_abs_value/2, max_abs_value),
    c("#2C75B3",      "#5B9BD5",       "#9DC3E6",       "#F9EAE1", "#F4B183",      "#ED7D31",      "#C55A11")
  )


  column_ha <- HeatmapAnnotation(
    Gene = c("TP53", "CTNNB1"),
    col = list(
      Gene = c("TP53" = "#D55E00", "CTNNB1" = "#0072B2")
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = col_fontsize)
  )


  ht <- Heatmap(
    heatmap_matrix,
    name = "Activity Score Change",
    col = col_fun,
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = col_fontsize, rot = 0),
    top_annotation = column_ha,
    cell_fun = function(j, i, x, y, width, height, fill) {
      p_val <- p_value_matrix[i, j]
      star_text <-
        if (p_val < 0.001) {
          "***"
        } else if (p_val < 0.01) {
          "**"
        } else if (p_val < 0.05) {
          "*"
        } else {
          ""
        }
      if (star_text != "") {
        grid.text(star_text, x, y, gp = gpar(fontsize = star_fontsize))
      }
    },
    column_title = "Effect of Mutant Genes on Agent States",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Activity Score Change",
      at = round(seq(-max_abs_value, max_abs_value, length.out = 7), 2),
      labels = round(seq(-max_abs_value, max_abs_value, length.out = 7), 2),
      legend_height = unit(4, "cm"),
      title_gp = gpar(fontsize = legend_title_fontsize),
      labels_gp = gpar(fontsize = legend_labels_fontsize)
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5)
  )


  print(ht)


  return(list(
    heatmap_data = heatmap_data,
    p_value_matrix = p_value_matrix
  ))
}


plot_pathway_comparison <- function(tp53_data, ctnnb1_data, pathways, output_dir, simplify_names = FALSE) {
  log_message("Plotting TP53/CTNNB1 mutation impact on specific pathways...
", print_to_console = TRUE)


  tp53_subset <- tryCatch({

    tp53_data %>%
      filter(AgentState %in% pathways) %>%
      mutate(
        Gene = "TP53",
        Change = Mean_Score_Mutant - Mean_Score_WT,
        Significance = ifelse(P_Value < 0.05, "Significant", "Non-significant")
      ) %>%
      select(AgentState, Gene, Change, P_Value, Significance)
  }, error = function(e) {
    log_message(paste0("dplyr filter/mutate/select failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


    result <- tp53_data[tp53_data$AgentState %in% pathways, ]
    result$Gene <- "TP53"
    result$Change <- result$Mean_Score_Mutant - result$Mean_Score_WT
    result$Significance <- ifelse(result$P_Value < 0.05, "Significant", "Non-significant")


    result <- result[, c("AgentState", "Gene", "Change", "P_Value", "Significance")]

    return(result)
  })

  ctnnb1_subset <- tryCatch({

    ctnnb1_data %>%
      filter(AgentState %in% pathways) %>%
      mutate(
        Gene = "CTNNB1",
        Change = Mean_Score_Mutant - Mean_Score_WT,
        Significance = ifelse(P_Value < 0.05, "Significant", "Non-significant")
      ) %>%
      select(AgentState, Gene, Change, P_Value, Significance)
  }, error = function(e) {
    log_message(paste0("dplyr filter/mutate/select failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


    result <- ctnnb1_data[ctnnb1_data$AgentState %in% pathways, ]
    result$Gene <- "CTNNB1"
    result$Change <- result$Mean_Score_Mutant - result$Mean_Score_WT
    result$Significance <- ifelse(result$P_Value < 0.05, "Significant", "Non-significant")


    result <- result[, c("AgentState", "Gene", "Change", "P_Value", "Significance")]

    return(result)
  })


  combined_data <- rbind(tp53_subset, ctnnb1_subset)


  if (simplify_names) {

    combined_data$AgentState <- gsub("KEGG_", "", combined_data$AgentState)
    combined_data$AgentState <- gsub("_PATHWAY", "", combined_data$AgentState)
    combined_data$AgentState <- gsub("_", " ", combined_data$AgentState)
  }


  pathway_order <- tryCatch({

    combined_data %>%
      group_by(AgentState) %>%
      summarise(mean_abs_change = mean(abs(Change))) %>%
      ungroup() %>%
      arrange(desc(mean_abs_change)) %>%
      pull(AgentState)
  }, error = function(e) {
    log_message(paste0("dplyr group_by/summarize/arrange/pull failed, using base R. Error: ", e$message, "\n"), print_to_console = TRUE)


    agg_data <- aggregate(abs(Change) ~ AgentState, data = combined_data, FUN = mean)
    colnames(agg_data)[2] <- "mean_abs_change"
    ordered_data <- agg_data[order(agg_data$mean_abs_change, decreasing = TRUE), ]
    return(ordered_data$AgentState)
  })

  combined_data$AgentState <- factor(combined_data$AgentState, levels = pathway_order)


  p <- ggplot(combined_data, aes(x = AgentState, y = Change, fill = Gene)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("TP53" = "#D55E00", "CTNNB1" = "#0072B2")) +
    geom_text(aes(label = ifelse(Significance == "Significant", "*", ""),
                  y = ifelse(Change >= 0, Change + 0.02, Change - 0.02)),
              position = position_dodge(width = 0.8),
              vjust = ifelse(combined_data$Change >= 0, 0, 1),
              size = 5) +
    coord_flip() +
    theme_publication(base_size = 14) +
    labs(
      title = "Comparison of TP53 and CTNNB1 Mutations on Metabolism Pathways",
      y = "Activity Score Change (Mutant - Wild Type)",
      x = "",
      fill = "Mutant Gene"
    )


  comparison_csv <- file.path(output_dir, "Figure6C_Metabolism_Pathway_Comparison_Data.csv")
  write.csv(combined_data, comparison_csv, row.names = FALSE)
  log_message(paste0("Pathway comparison data saved: ", comparison_csv, "\n"), print_to_console = TRUE)

  return(p)
}


plot_correlation <- function(tp53_data, ctnnb1_data, output_dir, highlight_pathways = NULL) {
  log_message("Plotting correlation scatter of TP53/CTNNB1 mutations on states...
", print_to_console = TRUE)


  tp53_fc <- tp53_data %>%
    mutate(Change_TP53 = Mean_Score_Mutant - Mean_Score_WT) %>%
    select(AgentState, Change_TP53, P_Value) %>%
    rename(P_Value_TP53 = P_Value)

  ctnnb1_fc <- ctnnb1_data %>%
    mutate(Change_CTNNB1 = Mean_Score_Mutant - Mean_Score_WT) %>%
    select(AgentState, Change_CTNNB1, P_Value) %>%
    rename(P_Value_CTNNB1 = P_Value)


  merged_data <- tp53_fc %>%
    inner_join(ctnnb1_fc, by = "AgentState") %>%
    mutate(
      Quadrant = case_when(
        P_Value_TP53 < 0.05 & P_Value_CTNNB1 < 0.05 & Change_TP53 > 0 & Change_CTNNB1 > 0 ~ "Both Up",
        P_Value_TP53 < 0.05 & P_Value_CTNNB1 < 0.05 & Change_TP53 < 0 & Change_CTNNB1 < 0 ~ "Both Down",
        P_Value_TP53 < 0.05 & P_Value_CTNNB1 < 0.05 & Change_TP53 > 0 & Change_CTNNB1 < 0 ~ "TP53 Up, CTNNB1 Down",
        P_Value_TP53 < 0.05 & P_Value_CTNNB1 < 0.05 & Change_TP53 < 0 & Change_CTNNB1 > 0 ~ "TP53 Down, CTNNB1 Up",
        P_Value_TP53 < 0.05 & Change_TP53 > 0 ~ "TP53 Up Only",
        P_Value_TP53 < 0.05 & Change_TP53 < 0 ~ "TP53 Down Only",
        P_Value_CTNNB1 < 0.05 & Change_CTNNB1 > 0 ~ "CTNNB1 Up Only",
        P_Value_CTNNB1 < 0.05 & Change_CTNNB1 < 0 ~ "CTNNB1 Down Only",
        TRUE ~ "Not Significant"
      ),
      Highlight = AgentState %in% highlight_pathways,
      Label = ifelse(AgentState %in% c("LipidMetabolism_Risk_LMRG", "KEGG_PROTEASOME"),
                     as.character(AgentState),
                     "")
    )


  corr <- cor.test(merged_data$Change_TP53, merged_data$Change_CTNNB1)
  lm_fit <- lm(Change_CTNNB1 ~ Change_TP53, data = merged_data)
  r_squared <- summary(lm_fit)$r.squared


  quadrant_colors <- c(
    "Both Up" = "#D55E00",
    "Both Down" = "#0072B2",
    "TP53 Up, CTNNB1 Down" = "#CC79A7",
    "TP53 Down, CTNNB1 Up" = "#009E73",
    "TP53 Up Only" = "#E69F00",
    "TP53 Down Only" = "#56B4E9",
    "CTNNB1 Up Only" = "#F0E442",
    "CTNNB1 Down Only" = "#999999",
    "Not Significant" = "grey80"
  )


  p <- ggplot(merged_data, aes(x = Change_TP53, y = Change_CTNNB1)) +

    geom_smooth(method = "lm", color = "grey50", linetype = "dashed", se = TRUE, alpha = 0.2) +

    geom_point(aes(color = Quadrant, size = -log10(pmin(P_Value_TP53, P_Value_CTNNB1))),
               alpha = ifelse(merged_data$Highlight, 1, 0.7)) +

    geom_point(data = subset(merged_data, Highlight),
              aes(color = Quadrant), size = 5, shape = 21, stroke = 2) +

    ggrepel::geom_text_repel(
      data = subset(merged_data, Label != ""),
      aes(label = Label),
      size = 5,
      box.padding = 0.8,
      segment.alpha = 0.8,
      force = 10,
      min.segment.length = 0,
      max.overlaps = Inf
    ) +

    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

    scale_color_manual(values = quadrant_colors) +
    scale_size_continuous(range = c(1, 5)) +

    labs(
      title = "Correlation of TP53 and CTNNB1 Mutation Effects",
      subtitle = sprintf("Pearson R = %.2f (P = %.2e), R² = %.2f",
                        corr$estimate, corr$p.value, r_squared),
      x = "Activity Change by TP53 Mutation (Log2FC)",
      y = "Activity Change by CTNNB1 Mutation (Log2FC)",
      color = "Change Pattern",
      size = "-log10(P-value)"
    ) +
    theme_publication(base_size = 22)


  correlation_csv <- file.path(output_dir, "Figure6F_Correlation_Data.csv")
  write.csv(merged_data, correlation_csv, row.names = FALSE)
  log_message(paste0("Correlation scatter data saved: ", correlation_csv, "\n"), print_to_console = TRUE)


  correlation_stats <- data.frame(
    correlation_coefficient = corr$estimate,
    p_value = corr$p.value,
    r_squared = r_squared,
    conf_int_lower = corr$conf.int[1],
    conf_int_upper = corr$conf.int[2]
  )
  stats_csv <- file.path(output_dir, "Figure6F_Correlation_Statistics.csv")
  write.csv(correlation_stats, stats_csv, row.names = FALSE)
  log_message(paste0("Correlation statistics saved: ", stats_csv, "\n"), print_to_console = TRUE)

  return(p)
}


main <- function(output_dir = "../results") {

  load_packages()


  data <- read_data()
  tp53_data <- data$tp53_data
  ctnnb1_data <- data$ctnnb1_data
  output_dir <- data$output_dir


  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_message(paste0("Creating output directory: ", output_dir, "\n"), print_to_console = TRUE)
  }


  tp53_significant <- filter_significant_results(tp53_data)
  ctnnb1_significant <- filter_significant_results(ctnnb1_data)


  write.csv(tp53_significant, file.path(output_dir, "Figure6_TP53_Significant_Results.csv"), row.names = FALSE)
  write.csv(ctnnb1_significant, file.path(output_dir, "Figure6_CTNNB1_Significant_Results.csv"), row.names = FALSE)
  log_message("Significant filtering result saved
", print_to_console = TRUE)


  results_dir <- output_dir
  log_message(paste0("Using output directory: ", results_dir, "\n"), print_to_console = TRUE)

  tryCatch({

    tp53_volcano <- plot_volcano(tp53_data, "TP53", output_dir)
    ctnnb1_volcano <- plot_volcano(ctnnb1_data, "CTNNB1", output_dir)


    ggsave(file.path(results_dir, "Figure6A_TP53_Volcano.pdf"), tp53_volcano, width = 10, height = 8)
    ggsave(file.path(results_dir, "Figure6B_CTNNB1_Volcano.pdf"), ctnnb1_volcano, width = 10, height = 8)
    log_message(paste0("Volcano plot saved to: ", results_dir, "\n"), print_to_console = TRUE)
  }, error = function(e) {
    log_message(paste0("Error plotting/saving volcano plot: ", e$message, "\n"), print_to_console = TRUE)
  })

  tryCatch({


    tp53_all_data <- tp53_significant %>%
      select(AgentState, Mean_Score_Mutant, Mean_Score_WT, P_Value) %>%
      rename(
        TP53_Mean_Score_Mutant = Mean_Score_Mutant,
        TP53_Mean_Score_WT = Mean_Score_WT,
        TP53_P_Value = P_Value
      )


    ctnnb1_all_data <- ctnnb1_significant %>%
      select(AgentState, Mean_Score_Mutant, Mean_Score_WT, P_Value) %>%
      rename(
        CTNNB1_Mean_Score_Mutant = Mean_Score_Mutant,
        CTNNB1_Mean_Score_WT = Mean_Score_WT,
        CTNNB1_P_Value = P_Value
      )


    all_pathways_data <- full_join(tp53_all_data, ctnnb1_all_data, by = "AgentState")


    write.csv(all_pathways_data, file.path(results_dir, "Figure6_All_Heatmap_Data.csv"), row.names = FALSE)
    log_message(paste0("All pathway heatmap raw data saved to:", results_dir, "/Figure6_All_Heatmap_Data.csv\n"),
                print_to_console = TRUE)


    write.csv(as.data.frame(heatmap_result$p_value_matrix), file.path(results_dir, "Figure6_All_Heatmap_PValues.csv"), row.names = TRUE)
    log_message(paste0("All pathway heatmap P-value matrix saved to:", results_dir, "/Figure6_All_Heatmap_PValues.csv\n"),
                print_to_console = TRUE)


    all_immune_pathways <- c(

      "T_cell_co_stimulation",
      "T_cell_co_inhibition",
      "APC_co_stimulation",
      "Check_point",
      "CCR",
      "Th1_cells",

      "CytokineSignaling_HighRisk_CSIRG"
    )


    tp53_all_immune <- tp53_data[tp53_data$AgentState %in% all_immune_pathways, ]
    ctnnb1_all_immune <- ctnnb1_data[ctnnb1_data$AgentState %in% all_immune_pathways, ]


    if (nrow(tp53_all_immune) > 0 || nrow(ctnnb1_all_immune) > 0) {
      log_message("Starting immune pathway heatmap...\n", print_to_console = TRUE)


      pdf(file.path(results_dir, "Figure6D_Key_Immune_Pathways_Heatmap.pdf"), width = 10, height = 8)
      heatmap_result <- plot_heatmap(
      tp53_all_immune, ctnnb1_all_immune, output_dir,
      top_n = length(all_immune_pathways),
      row_fontsize = 12,
      col_fontsize = 14,
      legend_title_fontsize = 12,
      legend_labels_fontsize = 10,
      star_fontsize = 12
    )
    dev.off()
      log_message("Immune pathway heatmap complete\n", print_to_console = TRUE)


      immune_raw_data_file <- file.path(results_dir, "Figure6D_Key_Immune_Pathways_Heatmap_Data.csv")
      write.csv(heatmap_result$heatmap_data, immune_raw_data_file, row.names = FALSE)
      log_message(paste0("Immune pathway heatmap raw data saved to:", immune_raw_data_file, "\n"),
                  print_to_console = TRUE)


      immune_pvalue_file <- file.path(results_dir, "Figure6D_Key_Immune_Pathways_Heatmap_PValues.csv")
      write.csv(as.data.frame(heatmap_result$p_value_matrix), immune_pvalue_file, row.names = TRUE)
      log_message(paste0("Immune pathway heatmap P-value matrix saved to:", immune_pvalue_file, "\n"),
                print_to_console = TRUE)


    } else {
      log_message("Warning: insufficient immune pathway data for heatmap.\n", print_to_console = TRUE)
    }
  }, error = function(e) {
    log_message(paste0("Error plotting/saving heatmap: ", e$message, "\n"), print_to_console = TRUE)
  })

  tryCatch({

    highlight_pathways <- c(
      "LipidMetabolism_Risk_LMRG",
      "KEGG_PROTEASOME"
    )


    correlation_plot <- plot_correlation(
      tp53_data = tp53_data,
      ctnnb1_data = ctnnb1_data,
      output_dir = output_dir,
      highlight_pathways = highlight_pathways
    )


    ggsave(
      file.path(results_dir, "Figure6F_Pathway_Correlation.pdf"),
      correlation_plot,
      width = 14,
      height = 10,
      dpi = 300
    )
    log_message(paste0("Correlation scatter with highlighted pathways saved to:", results_dir, "\n"),
                print_to_console = TRUE)
  }, error = function(e) {
    log_message(paste0("Error plotting/saving correlation scatter: ", e$message, "\n"),
                print_to_console = TRUE)
  })


  tryCatch({

    important_pathways <- c(
      "LipidMetabolism_Risk_LMRG",
      "Metabolism_General_Risk_MRG",
      "KEGG_FATTY_ACID_METABOLISM",
      "KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
      "KEGG_RETINOL_METABOLISM",
      "KEGG_TRYPTOPHAN_METABOLISM"
    )


    p_simplified <- plot_pathway_comparison(
      tp53_data = tp53_data,
      ctnnb1_data = ctnnb1_data,
      pathways = important_pathways,
      output_dir = output_dir,
      simplify_names = FALSE
    )


    simplified_pdf <- file.path(output_dir, "Figure6C_Metabolism_Pathway_Comparison.pdf")
    ggsave(simplified_pdf, p_simplified, width = 12, height = 8)
    log_message(paste0("Original-name pathway comparison plot saved to: ", simplified_pdf, "\n"), print_to_console = TRUE)

  }, error = function(e) {

    log_message(paste0("Error plotting/saving pathway comparison: ", e$message, "\n"), print_to_console = TRUE)
  })


  log_message(paste0("All saved violin plots were generated in Figure6.R, located at:", results_dir, "\n"), print_to_console = TRUE)
  log_message("All individual figures saved, skipping figure combination step\n", print_to_console = TRUE)
  log_message("All important-result visualizations complete!\n", print_to_console = TRUE)
}


main()
