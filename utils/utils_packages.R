
 
# CRAN和bioconductor
pkgs <- c('Seurat', 'tidyverse', 'signac', 'devtools', 'DropletUtils','',
          '','',''
          )

#devtools::install_local(package_name,force = T,quiet = F)
check_pkgs(pkgs)

#from https://github.com/compbioNJU/scPlant/blob/master/R/load_shinyApp.R  
check_pkgs <- function(pkgs) {
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  unavail <- c()
  for (pkg in pkgs) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      unavail <- c(unavail, pkg)
    }
  }
  if (length(unavail) > 0) {
    message("These packages are required but not installed: \n", paste(strwrap(paste(unavail, collapse = ", ")), collapse = "\n"))
      for (pkg in unavail) {
        tryCatch(install.packages(pkg),
                 error = function(e) {
                   message("Installation from CRAN failed: \n", e$message, "\nTry to install ", pkg, " from Bioconductor")
                   if(!requireNamespace("BiocManager", quietly = TRUE)) {
                     install.packages("BiocManager")
                   }
                   BiocManager::install(pkg)
                 })
    } else {
      stop("You need to manually install these packages.")
    }
  }
}




#github
github_list=c('samuel-marsh/scCustomize',
                'mojaveazure/seurat-disk',
                '' 
            )
