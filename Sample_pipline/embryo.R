project_name <- "Embryo"
pal <- viridis::viridis(n = 10)
source("/public/home/luoliheng/SINGLE/scRNA_analysis/utils.R", chdir = TRUE)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(SeuratWrappers))

work_dir <- "/public/home/luoliheng/embryonic_cells/new_analysis/trajectory/"

velo_data <- ReadVelocity(file = "/public/home/luoliheng/embryonic_cells/rawdata/velocyto/combined.loom") %>%
  as.Seurat()

velo_data$type <- colnames(velo_data) %>% str_extract("^[^:]+")

  velo_data<-  qc_check(velo_data,
    species = "hs",
    Find_doublet = FALSE
  )
  qc_plot%>% 
  

  SCTransform(assay = "spliced") %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20)
  
  
   %>%
  RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.02)

