project_name <- "MSC"
work_dir <- "/public/home/luoliheng/yujia/analysis/gz+jz/re_analysis_gz/"
pal <- viridis::viridis(n = 10)
source("/public/home/luoliheng/SINGLE/R/utils.R", chdir = TRUE)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(scCustomize))

all_data <- readRDS(save_file(fun = saveRDS, name_string = "qc"))

run_sctransform
require(clustree)
pdf(save_file(fun = pdf, name_string = "after_qc"))
gruop_var <- "type"
print(DimPlot(all_data, split.by = "orig.ident", ncol = 2))
print(DimPlot_scCustom(all_data, group.by = gruop_var, ggplot_default_colors = TRUE, figure_plot = TRUE))
print(DimPlot_scCustom(all_data, group.by = "Phase", ggplot_default_colors = TRUE, figure_plot = TRUE))
for (i in feats) {
    print(FeaturePlot_scCustom(all_data, colors_use = pal, features = i))
}
for (i in resolutions) {
    group <- paste0(assay_use, "_snn_res.", i)
    print(DimPlot(all_data, group.by = group, label = T))
}
print(clustree(all_data@meta.data, prefix = paste0(assay_use, "_snn_res.")))
dev.off()


require(readxl)

Idents(all_data) <- paste0(assay_use, "_snn_res.", 0.2)
markers <- read_excel("/public/home/luoliheng/yujia/analysis/gz+jz/down_stream_gz/files/markers_idents.xlsx", sheet = 1)

pdf(save_file(fun = pdf, name_string = "markers"))

plot_makers(unfilterd_markers = markers, all_data, dot_plot = TRUE, max_markers = 40, cluster = FALSE, feature_plot = TRUE, col_dot = viridis_plasma_dark_high, col_fea = pal)

dev.off()
