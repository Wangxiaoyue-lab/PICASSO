work_dir <- "/public/home/luoliheng/yujia/analysis/gz+jz/re_analysis_jz/"
project_name <- "MSC"
source("/public/home/luoliheng/SINGLE/scRNA_analysis/utils.R", chdir = TRUE)

pal <- viridis::viridis(n = 10)

namelist <- read.table("/public/home/luoliheng/yujia/counts_results/namelist", header = F)

all_list <- namelist[13:18, ] %>%
    map(~ {
        paste0("/public/home/luoliheng/yujia/counts_results/count_", .x, "/outs/filtered_feature_bc_matrix") %>%
            Read10X(data.dir = .) %>%
            CreateSeuratObject(min.cells = 10, min.features = 200, project = .x)
    })
all_list[1:3] <- lapply(all_list[1:3], function(x) {
    x$type <- "NC"
    return(x)
})

all_list[4:6] <- lapply(all_list[4:6], function(x) {
    x$type <- "HF"
    return(x)
})

save_file(fun = pdf, name_string = "before_qc")
all_data <- merge(all_list[[1]], all_list[-1])
all_data <- all_data %>%
    qc_check(
        species = "mm",
        Find_doublet = TRUE
    ) %>%
    qc_check_plot(feature_scatter = TRUE, vln_group = "orig.ident", dim_group = c("orig.ident", "type", "Phase"), feats = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb"), colors = pal)

dev.off()

save_file(data = all_data, fun = saveRDS, name_string = "qc")


run_sctransform <- TRUE
all_data <- readRDS(save_file(fun = saveRDS, name_string = "qc"))
all_data <- qc_process(all_data,
    dim_use = 20,
    resolutions = c(0.1, 0.2, 0.3, 0.5),
    run_harmony = TRUE,
    max_mt = 20,
    max_hb = 5,
    nFeature_RNA_min = 500,
    nFeature_RNA_max = 7500,
    vars_to_regress = c("percent_mito", "S.Score", "G2M.Score"),
    group_by_vars = "orig.ident"
)

save_file(data = all_data, fun = saveRDS, name_string = "qc")
