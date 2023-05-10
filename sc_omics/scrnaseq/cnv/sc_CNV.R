sc_cnv <- function(object, method, output) {

}


sc_cnv_infercnv <- function(object,
                            outpath,
                            celltype_col, # metadta里标记细胞类型的列的列名
                            normal_cell, # 作为参考的正常细胞
                            type, # 按样本预测cnv还是按亚克隆预测cnv
                            denoise = T,
                            HMM = T,
                            ncores = NULL) {
    # 0 check the settings
    library(infercnv)
    library(ggplot2)
    assertthat::assert_that(class(object) == "Seurat")
    species <- check_species(object)
    assertthat::assert_that(species %in% c("human", "mouse"))
    assertthat::assert_that(type %in% c("samples", "subclusters"))
    assertthat::assert_that(all(celltype_col %in% colnames(object@meta.data)))
    assertthat::assert_that(all(normal_cell %in% unique(object@meta.data[, celltype_col])))
    ncores <- ncores %||% 20
    # 1 prepare the input
    dir.create(paste0(outpath, "/input"))
    dir.create(paste0(outpath, "/output"))
    dir.create(paste0(outpath, "/output", "/picture"))
    dir.create(paste0(outpath, "/output", type))
    gene2cnv <- rownames(object) %>%
        AnnoProbe::annoGene(., "SYMBOL", species) %>%
        arrange(chr, start) %>%
        select(SYMBOL, chr, start, end) %>%
        distinct(SYMBOL, .keep_all = T)
    write.table(gene2cnv,
        file = paste0(outpath, "/input", "/gene_file.txt"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = F
    )
    exp2cnv <- FetchData(object,
        vars = gene2cnv$SYMBOL,
        cells = colnames(object),
        layer = "count"
    )
    write.table(exp2cnv,
        file = paste0(outpath, "/input", "/exp_file.txt"),
        sep = "\t",
        quote = F
    )
    anno2cnv <- object@meta.data %>%
        mutate(cellnames2cnv = rownames(.)) %>%
        select(cellnames2cnv, !!sym(celltype_col))
    write.table(anno2cnv,
        file = paste0(outpath, "/input", "/anno_file.txt"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = F
    )
    rm(list = c("exp2cnv", "gene2cnv", "anno2cnv"))
    message("Succeed to prepare the input")
    # 2 run the program
    future::plan("multiprocess", workers = ncores)
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = paste0(outpath, "/input", "/exp_file.txt"), # raw counts matrix
        annotations_file = paste0(outpath, "/input", "/anno_file.txt"), # an annotations file which indicates which cells are tumor vs. normal.
        delim = "\t",
        gene_order_file = paste0(outpath, "/input", "/gene_file.txt"), # gene/chromosome positions file
        ref_group_names = normal_cell # reference cell name
    )
    if (type == "samples") {
        infercnv_obj <- infercnv::run(infercnv_obj,
            cutoff = 0.1,
            # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
            out_dir = paste0(outpath, "/output", type),
            cluster_by_groups = TRUE,
            analysis_mode = type,
            denoise = TRUE,
            HMM = TRUE
        )
    } else {
        infercnv_obj <- infercnv::run(infercnv_obj,
            cutoff = 0.1,
            out_dir = paste0(outpath, "/output", type),
            cluster_by_groups = F,
            analysis_mode = type,
            hclust_method = "ward.D2",
            tumor_subcluster_partition_method = "random_trees",
            tumor_subcluster_pval = 0.05,
            denoise = TRUE,
            HMM = TRUE
        )
    }
    saveRDS(infercnv_obj, file = paste0(outpath, "/output", "/infercnv_obj.RDS"))
    infercnv::plot_cnv(infercnv_obj,
        out_dir = outpath,
        output_filename = "infercnv.plot",
        plot_chr_scale = F, # 是否画染色体全长
        output_format = "pdf",
        x.range = "auto",
        # x.center = 1,
        title = "infercnv",
        color_safe_pal = FALSE
    )
    return(infercnv_obj)
}

sc_cnv_copykat <- function(object,
                           outpath,
                           celltype_col, # metadta里标记细胞类型的列的列名
                           normal_cell, # 作为参考的正常细胞
                           ncores) {
    # 0 check the setting
    library(copykat)
    assertthat::assert_that(all(celltype_col %in% colnames(object@meta.data)))
    assertthat::assert_that(all(normal_cell %in% unique(object@meta.data[, celltype_col])))
    ncores <- ncores %||% 20
    # 1 prepare the input
    dir.create(paste0(outpath, "/input"))
    dir.create(paste0(outpath, "/output"))
    dir.create(paste0(outpath, "/output", "/picture"))
    exp2cnv <- FetchData(object,
        vars = gene2cnv$SYMBOL,
        cells = colnames(object),
        layer = "count"
    )
    write.table(exp2cnv,
        file = paste0(outpath, "/input", "/exp_file.txt"),
        sep = "\t",
        quote = F,
        row.names = T,
        col.names = T
    )
    anno2cnv <- object@meta.data %>%
        mutate(cellnames2cnv = rownames(.)) %>%
        dplyr::select(cellnames2cnv, !!sym(celltype_col))
    normal2cnv <- anno2cnv %>% 
        dplyr::filter(!!sym(celltype_col) %in% normal_cell) %>%
        dplyr::pull(cellnames2cnv)
    assertthat::assert_that(length(normal2cnv) > 50)
    # 2 run the program
    copykat_obj <- copykat(
        rawmat = exp2cnv,
        id.type = "S", # id是gene symbol
        cell.line = "no",
        ngene.chr = 5,
        win.size = 25,
        KS.cut = 0.1,
        norm.cell.names = normal2cnv,
        sam.name = "copykat_", # 文件前缀
        distance = "euclidean",
        n.cores = ncores
    )
    saveRDS(copykat_obj, file = paste0(outpath, "/output", "/copykat_obj.RDS"))
    # 肿瘤是否恶性
    pred.test <- data.frame(copykat.test$prediction)
    pred.test <- lefj_join(x = pred.test, y = anno2cnv, 
        dplyr::join_by(cell.names == cellnames2cnv)) #copykat会自动过滤基因数少于100的细胞
    saveRDS(pred.test, file = paste0(outpath, "/output", "/pred.test.RDS"))
    print(table(pred.test[, c(celltype_col, "copykat.pred")]))
    # cnv事件坐标与检测量
    CNA.test <- data.frame(copykat.test$CNAmat)
    saveRDS(CNA.test, file = paste0(outpath, "/output", "/CNA.test.RDS"))
    # 3 plot the result
    ## chromosome color
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
    chr <- as.numeric(CNA.test$chrom) %% 2+1
    rbPal1 <- colorRampPalette(c('black','grey'))
    CHR <- rbPal1(2)[as.numeric(chr)]
    chr1 <- cbind(CHR,CHR)
    ## scale color
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    com.preN <- pred.test$copykat.pred
    pred <- rbPal5(2)[as.numeric(factor(com.preN))]
    cells <- rbind(pred,pred)
    col_breaks = c(seq(-1,-0.4,length=50),
                   seq(-0.4,-0.2,length=150),
                   seq(-0.2,0.2,length=600),
                   seq(0.2,0.4,length=150),
                   seq(0.4, 1,length=50))
    pdf(paste0(outpath, "/output", "/picture", "/CNA_tumor_normal.pdf"),
        width = 14, height = 16)
    heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),
            dendrogram="r", 
            distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"),
            hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,
            RowSideColors=cells,
            Colv=NA, Rowv=TRUE,
            notecol="black",
            col=my_palette,
            breaks=col_breaks, 
            key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
    legend("topright", 
        paste("pred.",names(table(com.preN)),sep=""), 
        pch=15,
        col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], 
        cex=0.6, bty="n")
    dev.off()
    # 肿瘤根据层次聚类分亚群
    tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
    tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
    hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), 
        method = "ward.D2")
    hc.umap <- cutree(hcc,2)
    rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
    subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
    cells <- rbind(subpop,subpop)
    heatmap.3(t(tumor.mat),dendrogram="r", 
        distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), 
        hclustfun = function(x) hclust(x, method="ward.D2"),
        ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
        notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
        keysize=1, density.info="none", trace="none",
        cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
        symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')

    return(copykat.test)
}


sc_cnv_score <- function() {
    # https://zhuanlan.zhihu.com/p/433234064
    # https://www.jianshu.com/p/1fa1fd4f97ff
}







sc_cnv_honeybadger <- function() {
    library(HoneyBADGER)
    # https://jef.works/HoneyBADGER/
    ## 表现不太好，使用不多，有需要再补
}

sc_cnv_casper <- function() {
    library(CaSpER)
    # https://github.com/akdess/CaSpER/blob/master/demo/MM135_10X.R
}
