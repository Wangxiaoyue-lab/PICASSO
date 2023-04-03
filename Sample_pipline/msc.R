work_dir <- "/public/home/luoliheng/yujia/analysis/gz+jz/re_analysis_gz/markers/"
project_name <- "HSC"
source("/public/home/luoliheng/SINGLE/scRNA_analysis/utils.R", chdir = TRUE)
pal <- viridis::viridis(n = 10)

file_dir <- "/public/home/luoliheng/yujia/analysis/gz+jz/marker_files/marker_0330.xlsx"

require(readxl)
marker <- read_excel(file_dir, sheet = 1)
all_data <- readRDS("/public/home/luoliheng/yujia/analysis/gz+jz/re_analysis_gz/MSC_qc.rds")
Idents(all_data) <- "SCT_snn_res.0.2"
save_file(fun = pdf, name_string = "markers0330")
plot_makers(marker, all_data, dot_plot = TRUE, max_markers = 40, cluster = "no", feature_plot = FALSE, col_dot = viridis_plasma_dark_high, col_fea = pal)
dev.off()
