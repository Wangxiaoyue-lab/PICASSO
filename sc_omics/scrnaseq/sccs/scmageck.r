#######################################
# DATE； 20230320
# AUTHOR: Cao Jun
# DESCRIPTION:运行scMAGeCK的RRA和LR
######################################
# 加载包

library(scMAGeCK)



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
## 是否运行scmageckLR
runlr <- T



## 靶点选择
load(paste0(dir, "/data/targets"))
targets <- corper_fen_anno %>%
    # filter(perturb_kmeans %in% c(2,4)) %>% pull(gene)
    filter(perturb_kmeans %in% c(1, 2, 3, 4)) %>% pull(gene)

## 表达数据
if (T) {
    ### 只输入需要分析的扰动的细胞
    data.seurat <- readRDS(file = paste0(dir, exp))
    data.DM <- data.seurat %>%
        subset(subset = type == "DM") %>%
        .[, .$gene %in% c(targets, "Non-Targeting")] %>%
        NormalizeData() %>%
        FindVariableFeatures(nfeatures = 4000) %>%
        ScaleData()
    saveRDS(data.DM, file = paste0(dir, output, "/my_dm_filt.RDS"))
}

## metadata重命名cellnames,barcode,sgrna,gene,read_count,umi_count

# data.seurat=readRDS(file=paste0(dir,output,'/my_dm_filt.RDS'))
data.seurat <- readRDS(file = paste0(dir, "/data/DM_seurat_rds")) # 软链接


#################################################################
# scMAGeCK
## 将高变基因4000都进行scMAGeCK
if (runlr) {
    ## 1 构造元数据
    ### 元数据需要细胞名、sgrna名、sgrna序列、靶基因名、read计数和umi计数
    BARCODE <- data.seurat@meta.data %>%
        rownames_to_column("cellnames") %>%
        select(cellnames, barcode, sgrna, gene, read_count, umi_count) %>%
        mutate(read_count = 1, umi_count = 1, gene = str_replace(.$gene, "Non-Targeting", "non.target")) %>%
        mutate(barcode = ifelse(grepl("non_target", barcode), "non.target", barcode)) %>%
        rename(cell = cellnames)
    write.table(BARCODE, file = paste0(dir, output, "/BARCODE_sub_for_scmageck.txt"), quote = F, sep = "\t", col.names = T)


    ## 2 Run scMAGeCK lr

    RDS_path <- paste0(dir, output, "/my_dm_filt.RDS")
    BARCODE_path <- paste0(dir, output, "/BARCODE_sub_for_scmageck.txt")
    RRAPATH <- NULL
    NEGCTRL <- "non.target"

    lr_result <- scmageck_lr(BARCODE = BARCODE_path, RDS = RDS_path, LABEL = "scmageck_lr", NEGCTRL = NEGCTRL, PERMUTATION = 1000, SAVEPATH = NULL, LAMBDA = 0.01, GENE_FRAC = 0)
    lr_score <- lr_result[1][[1]]
    lr_score_pval <- lr_result[2][[1]]
    # save(lr_result,lr_score,lr_score_pval,file=paste0(dir,output,"/scmageck_lr.Rdata"))
    save(lr_result, lr_score, lr_score_pval, file = paste0(dir, output, "/p1234/scmageck_lr.Rdata"))
}

## 3 Run scMaGeCK RRA
# if(F){
# data.DM <- readRDS(file=paste0(dir,'/oscar_dm_filt.RDS'))
# scmageck_all=function(x){
# result <- scmageck_rra(GENE=x,BARCODE=BARCODE_path, RDS=RDS_path, LABEL='scmageck_rra', NEGCTRL = NEGCTRL, SAVEPATH=NULL, KEEPTMP=F,PATHWAY=F)
# return(result)
# }
# rra_result=lapply(VariableFeatures(data.DM),scmageck_all)
# names(rra_result) <- VariableFeatures(data.DM)
# s.DM <- readRDS(file=paste0(dir,'/oscar_dm_filt.RDS'))
# scmageck_all=function(x){
# result <- scmageck_rra(GENE=x,BARCODE=BARCODE_path, RDS=RDS_path, LABEL='scmageck_rra', NEGCTRL = NEGCTRL, SAVEPATH=NULL, KEEPTMP=F,PATHWAY=F)
# return(result)
# }
# rra_result=lapply(VariableFeatures(data.DM),scmageck_all)
# names(rra_result) <- VariableFeatures(data.DM)
# save(rra_result,file=paste0(dir,"/scmageck_rra_oscar.Rdata"))
# }
