source("/public/home/luoliheng/SINGLE/R/utils.R", chdir = TRUE)

suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(Seurat))
work_dir <- "/public/home/luoliheng/yujia/analysis/gan/basic-analysis/qc/"
project_name <- "LIVER"
species <- "mm"
all_data <- readRDS(paste0(work_dir, "gan_qc.rds"))

marker_genes <- find_markers(all_data, all = FALSE, ident = NULL, loop_var = 0:3, ident.1 = "HF", group.by = "type", only.pos = FALSE)

save_file(data = marker_genes, fun = write.csv, name_string = "markers_sz", row.names = FALSE, col.names = TRUE)
