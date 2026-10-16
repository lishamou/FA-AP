
library(ggplot2)
library(RColorBrewer)

sci_theme_red <- function() {
  theme_minimal(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, colour = "black"),
      legend.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 14, lineheight = 1.2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#BDBDBD", linewidth = 0.2),
      panel.border = element_rect(colour = "#BDBDBD", fill = NA, linewidth = 0.5),
      panel.background = element_blank(),
      plot.margin = margin(1, 2, 1, 1, "cm")
    )
}


sci_red_palette <- colorRampPalette(c("#67001f", "#E41A1C", "#f7f7f7", "#F8766D", "#67001f"))
