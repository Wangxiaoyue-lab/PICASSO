#######################################
# DATE； 20230320
# AUTHOR: Cao Jun
# DESCRIPTION:运行MIMOSCA
######################################
# 加载包

library(Seurat)
library(dplyr)
library(purrr)
library(stringr)
library(conflicted)
library(reticulate)


source("/public/home/caojun/project/OSCAR/0_trans_data/mimosca/R_mimosca_packages.R")
source("/public/home/caojun/project/OSCAR/0_trans_data/mimosca/R_mimosca_support.R")
conflict_scout()
conflict_prefer("filter", "dplyr")

# 调试
# srun -p normal --pty /public/home/caojun/anaconda3/envs/r4.2/bin/R

############################################################
# 加载变量
## task路径
dir <- "/public/home/caojun/project/OSCAR/small_bench_20230316"
## seurat对象，rds
exp <- "/data/DM_EM_seurat_rds"
### dm+em，没有normalization
## output 输出文件夹
output <- "/output/my/"
## 是否运行
run_mimosca <- F
run_newcoef <- T

## 靶点选择
load(paste0(dir, "/data/targets"))
targets <- corper_fen_anno %>%
    # dplyr::filter(perturb_kmeans %in% c(2,4)) %>% pull(gene)
    filter(perturb_kmeans %in% c(1, 2, 3, 4)) %>% pull(gene)


## 表达数据
if (F) {
    data.seurat <- readRDS(file = paste0(dir, exp))
    data.DM <- data.seurat %>%
        subset(subset = type == "DM") %>%
        .[, .$gene %in% c(targets, "Non-Targeting")] %>%
        FindVariableFeatures(nfeatures = 4000) %>%
        ScaleData()
    saveRDS(data.DM, file = paste0(dir, "/oscar_dm_filt.RDS"))
    # 未normalization！！！
}
# data.seurat <- readRDS(   )

###########################################################


# 运行MIMOSCA的R语言版本，获得扰动概率
if (run_mimosca) {
    ## 制作Y矩阵

    ## 制作X矩阵

    ## 嵌入sgrna信息

    ## 进行扰动概率推断
    # ko_inf函数
}

########################################################

# 截断后重新分配身份并获得新系数
if (run_newcoef) {
    # 从我本来运行结果里获得已经运行完的部分
    load(paste0(dir, "/data/seurat_and_prob"))
    data.seurat <- data.DM %>% NormalizeData()

    # load(  其他数据   ko_inf + data.seurat)

    # ko_inf是函数infer_ko_probs的输出
    ## 包括了models、coefficients、squared_error、ko_probabilities四个部分
    # 第一次运行的系数是不能用的，需要根据概率截断后重新算一遍系数
    binarize <- function(x, threshold = NA) {
        matd <- x
        if (is.na(threshold)) {
            threshold <- min(x) + (max(x) - min(x)) / 2
            print(paste("Threshold: ", threshold))
        }
        matd[matd <= threshold] <- 0
        matd[matd > threshold] <- 1
        matd
    }


    sgrna_prob <- ko_inf[["ko_probabilities"]] %>% binarize(threshold = 0.2)
    genes_target <- colnames(sgrna_prob) %>%
        stringr::str_split(pattern = "_", simplify = T) %>%
        .[, 1] %>%
        unique()
    gene_prob <- lapply(genes_target, function(x) {
        target_guide <- sgrna_prob %>%
            as.matrix() %>%
            as_tibble() %>%
            dplyr::select(contains(paste0(x, "_"))) %>%
            rowSums()
    }) %>%
        do.call(cbind, .) %>%
        as.data.frame()
    colnames(gene_prob) <- genes_target
    rownames(gene_prob) <- rownames(sgrna_prob)

    if (F) {
        # 比较以0.2为界相比我之前策略的差别
        ## 0.2的细胞会更多
        gene_prob %>%
            rowSums() %>%
            table()
        ## 0     1
        ## 4183 26034

        ## HMM的细胞会少一些
        cell_P %>%
            do.call(c, .) %>%
            length()
        # [1] 17406
    }

    # 重新计算P2P4基因的扰动系数
    if (T) {
        ## 制作X矩阵
        X <- gene_prob %>%
            dplyr::select(all_of(targets)) %>%
            dplyr::filter(rowSums(.) > 0)


        ## 制作Y矩阵 4000基因
        top4000 <- data.seurat %>%
            FindVariableFeatures(nfeatures = 4000) %>%
            VariableFeatures()
        Y <- data.seurat[["RNA"]]@data[top4000, rownames(X)] %>% t()

        ## 制作C矩阵
        C <- data.seurat@meta.data[rownames(X), c("nFeature_RNA", "Phase")]

        ## 进行多元线性回归
        XC <- reformulate(termlabels = c(colnames(C), colnames(X))) %>%
            model.matrix(., cbind(C, X)) %>%
            .[, -1]
        fits <- fit_glmnet(XC, Y, alpha = 0.5, family = "gaussian")
        coefs <- fits %>%
            map_par(get_coef, verbose = T) %>%
            do.call(cbind, .) %>%
            as.data.frame()
        colnames(coefs) <- names(fits)
        # save(coefs,file=paste0(dir,output,'/mimosca_4000.Rdata'))
        save(coefs, file = paste0(dir, output, "/p1234/mimosca_4000.Rdata"))
    }
}
