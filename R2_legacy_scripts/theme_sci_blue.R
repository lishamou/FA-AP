
library(ggplot2)
library(RColorBrewer)

sci_theme_blue <- function() {
  theme_minimal(base_family = "sans") +
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


sci_blue_palette <- colorRampPalette(c("#053061", "#377EB8", "#f7f7f7", "#56B1F7", "#08306b"))
