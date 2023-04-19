#######################################
# DATE； 20230321
# AUTHOR: Cao Jun
# DESCRIPTION:OSCAR策略
######################################

# options(digits=10)
# 加载包
library(Seurat)
library(dplyr)
library(purrr)
library(stringr)
library(SCENIC)
library(SCopeLoomR)
library(conflicted)
source("/public/home/caojun/project/OSCAR/small_bench_20230316/script/1_my_data_compare3/design4/utils.r")
# conflict_scout()
conflict_prefer("filter", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("slice", "dplyr")
# 调试
# srun -p normal --pty /public/home/caojun/anaconda3/envs/r4.2/bin/R



# 加载变量
## task路径
dir <- "/public/home/caojun/project/OSCAR/small_bench_20230316"
## seurat对象，rds
exp <- "/data/DM_EM_seurat_rds"
### dm+em，没有normalization
## output 输出文件夹
output <- "/output/my/"
## 是否运行
run_pySCENIC <- F
# 物种
species <- "mm" #' hs'

effect_size <- T




## 靶点选择
load(paste0(dir, "/data/targets"))
targets <- corper_fen_anno %>%
    # filter(perturb_kmeans %in% c(2,4)) %>% pull(gene)
    filter(perturb_kmeans %in% c(1, 2, 3, 4)) %>% pull(gene)
## 表达数据
data.seurat <- readRDS(file = paste0(dir, exp)) %>% NormalizeData()
## dm+em，没有normalization


#################################################################
## 运行pySCENIC
## 先dm和em一起运行
if (run_pySCENIC) {
    # R输出csv文件
    data.seurat[["RNA"]]@counts %>%
        as.matrix() %>%
        t() %>%
        write.csv(file = paste0(dir, output, "count_for_pyscenic.csv"))

    # python转写为loom文件
    use_python("/public/home/caojun/anaconda3/envs/mamba_py3.8/bin/python")
    py_run_file(paste0(dir, "/script/csv2loom.py"), local = list(dir = dir, project = output, output = output))

    # pyscenic运行
    arg1 <- dir
    arg2 <- ouput # input_loom
    arg3 <- species # species: hs/mm
    paste("sh pyscenic.sh", arg1, arg2, arg3) %>% system()
}

#####################################################################

# 整理scenic结果获得OSCAR矩阵

if (effect_size) {
    ## regulon数据
    loom_path <- paste0(dir, "/data/my_scenic")

    regulonAUC <- open_loom(loom_path) %>%
        get_regulons_AUC(column.attr.name = "RegulonsAUC")
    regulons <- open_loom(loom_path) %>%
        get_regulons(column.attr.name = "Regulons") %>%
        regulonsToGeneLists()
    AUC_value <- regulonAUC@assays@data@listData$AUC

    ## oscar效应矩阵
    meta_OSCAR <- data.seurat@meta.data %>%
        select(barcode, gene, type) %>%
        filter(type == "DM", gene %in% c(targets, "Non-Targeting")) %>%
        rename(sgrna = barcode)

    mean_target <- function(meta, AUC, target) {
        # @output:The regulons X Perturbations matrix
        # @meta:The meta.data of Seurat object (includes gene and sgrna)
        # @AUC:The matrix of AUC value
        # @target:The choice to calulate the mean(gene or sgRNA)
        AUC <- AUC %>%
            as.data.frame() %>%
            select(all_of(rownames(meta)))
        pass_sgrna <- meta %>%
            group_by(sgrna) %>%
            summarize(count = n()) %>%
            filter(count >= 40) %>%
            pull(sgrna)
        meta <- meta %>%
            filter(sgrna %in% pass_sgrna) %>%
            rename_with(~ gsub(target, "group", .x), starts_with(target))
        mean_regulon <- lapply(meta$group %>% unique(), function(x) {
            cell_x <- meta %>%
                filter(group == x) %>%
                rownames()
            AUC_x <- AUC %>%
                select(all_of(cell_x)) %>%
                map_df(as.numeric) %>%
                rowMeans()
        }) %>%
            do.call(cbind, .) %>%
            as_tibble() %>%
            rename_all(~ meta$group %>% unique()) %>%
            mutate(Genes = rownames(AUC)) %>%
            select(Genes, everything())

        return(mean_regulon)
    }

    mean_target_OSCAR <- mean_target(meta = meta_OSCAR, AUC = AUC_value, target = "gene") %>% select(Genes, `Non-Targeting`, everything())

    effect_matrix <- mean_target_OSCAR %>%
        mutate(across(3:ncol(.), ~ . - mean_target_OSCAR[[2]])) %>%
        select(-2) %>%
        as.data.frame()

    effect_matrix_scale <- mean_target_OSCAR %>%
        mutate(across(3:ncol(.), ~ . - mean_target_OSCAR[[2]])) %>%
        select(-2) %>%
        mutate(across(2:ncol(.), ~ (. - mean(.)) / sd(.))) %>%
        as.data.frame()
    # save(effect_matrix,effect_matrix_scale,file=paste0(dir,output,"/mydata_oscar_effectsize.Rdata"))
    save(effect_matrix, effect_matrix_scale, file = paste0(dir, output, "/p1234/mydata_oscar_effectsize.Rdata"))
}
