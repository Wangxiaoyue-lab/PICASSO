```r
pkgs <- c('shiny', 'shinydashboard', 'data.table', 'DT', 'ggplot2', 'shinycssloaders', 'shinyWidgets',
            'dplyr', 'plyr', 'magrittr', 'Seurat', 'RColorBrewer', 'pheatmap', 'igraph', 'ggraph', 'circlize', 'ggpubr', 'cowplot',
            'factoextra', 'scatterpie',  'grid', 'topicmodels', 'tidytext',
            'ggwordcloud', 'wordcloud', 'pals', 'networkD3', 'echarts4r', 'ggnewscale')
  check_pkgs(pkgs)
check_pkgs <- function(pkgs) {
  unavail <- c()
  for (pkg in pkgs) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      unavail <- c(unavail, pkg)
    }
  }
  if (length(unavail) > 0) {
    if(!interactive()) {
      stop("These packages are required but not installed, you need to manually install: \n",
           paste(strwrap(paste(unavail, collapse = ", ")), collapse = "\n"))
    }
    message("These packages are required but not installed: \n", paste(strwrap(paste(unavail, collapse = ", ")), collapse = "\n"))
    answer = readline("Do you want to install them? [y|n] ")

    if(tolower(answer) %in% c("y", "yes")) {
      for (pkg in unavail) {
        tryCatch(install.packages(pkg),
                 error = function(e) {
                   message("Installation from CRAN failed: \n", e$message, "\nTry to install ", pkg, " from Bioconductor")
                   if(!requireNamespace("BiocManager", quietly = TRUE)) {
                     install.packages("BiocManager")
                   }
                   BiocManager::install(pkg)
                 })
      }
    } else {
      stop("You need to manually install these packages.")
    }
  }
}
```
