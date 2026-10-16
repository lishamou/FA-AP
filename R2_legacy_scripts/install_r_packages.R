

packages <- c(
  "tidyverse",
  "ggplot2",
  "reshape2",
  "gridExtra",
  "scales",
  "viridis",
  "fmsb",
  "showtext",
  "Cairo"
)


new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]


if(length(new_packages)) {
  cat("Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "https://cloud.r-project.org")
} else {
  cat("All required packages already installed.
")
}


cat("\nVerifying package loading:
")
for(pkg in packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat(sprintf("✓ %s  loaded successfully
", pkg))
  }, error = function(e) {
    cat(sprintf("✗ %s  failed to load: %s\n", pkg, e$message))
  })
}
