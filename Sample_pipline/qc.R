# Load in the utils.R script to use the functions defined in it
source("/public/home/luoliheng/SINGLE/R/utils.R", chdir = TRUE)

# Load in the tidyverse package and suppress the package startup messages
suppressPackageStartupMessages(require(tidyverse))
# Load in the Seurat package and suppress the package startup messages
suppressPackageStartupMessages(require(Seurat))

# Define the working directory where the data is stored
work_dir <- "/public/home/luoliheng/yujia/analysis/gan/basic-analysis/qc/"
# Define the project name
project_name <- "LIVER"
# Define the species
species <- "mm"
# Read in the data
all_data <- readRDS(paste0(work_dir, "gan_qc.rds"))

# Find the markers for each cluster
log_info("Finding markers...")
marker_genes <- find_markers(all_data, all = FALSE, ident = NULL, loop_var = 0:3, ident.1 = "HF", group.by = "type", only.pos = FALSE)
log_info("Finding markers done")

# Save the markers to a file
log_info("Saving markers to file...")
save_file(data = marker_genes, fun = write.csv, name_string = "markers_sz", row.names = FALSE, col.names = TRUE)
log_info("Saving markers to file done")
