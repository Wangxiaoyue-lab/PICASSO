sc_cnv_function <- function(...) {
    cat("# function list\n")
    cat("#-- sc_cnv_infercnv: Detection of CNV by infercnv\n")
    cat("#-- sc_cnv_copykat:  Detection of CNV by copykat\n")
    cat("#-- sc_cnv_honeybadger:  Detection of CNV by honeybadger\n")
    cat("#-- sc_cnv_casper:  Detection of CNV by casper\n")
    cat("#-- sc_cnv_score: Scoring of cells' CNV\n")
}


#' Detection of CNV in scRNA-seq Data
#'
#' This function uses the inferCNV package to detect Copy Number Variations (CNVs) in Single-Cell RNA sequencing data.
#'
#' @param object A 'Seurat' object.
#' @param task_path Output directory for saving results.
#' @param celltype_col The name of the column in metadata identifying cell types.
#' @param normal_cell The name of the normal cell type.
#' @param type The level predicting tbhe cnv(samples or subclusters).
#' @param denoise Use the denoise method for reducing noise in the data (default is TRUE).
#' @param HMM Use the Hidden Markov Model for smoothing the data (default is FALSE).
#' @param ncores The number of cores used in parallel (default is 20).
#'
#' @return A matrix of CNV scores corresponding to the input Seurat object.
#'
#' @examples
#' library(Seurat)
#' data("pbmc_small")
#' cnv_scores <- sc_cnv_infercnv(
#'     seurat_obj = object, task_path = "cnv_results",
#'     celltype_col = "seurat_clusters", normal_cell = "CD8T",
#'     denoise = TRUE, HMM = FALSE
#' )
#'
#' @import Seurat
#' @import inferCNV
#' @export
sc_cnv_infercnv <- function(object,
                            task_path,
                            prefix,
                            celltype_col,
                            normal_cell,
                            type,
                            denoise = T,
                            HMM = T,
                            ncores = NULL,
                            generate = F) {
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
    if (!dir.exists(task_path)) {
        dir.create(task_path)
    }
    if (!dir.exists(paste0(task_path, "/input"))) {
        dir.create(paste0(task_path, "/input"))
    }
    if (!dir.exists(paste0(task_path, "/output"))) {
        dir.create(paste0(task_path, "/output"))
    }
    if (!dir.exists(paste0(task_path, "/output", "/picture"))) {
        dir.create(paste0(task_path, "/output", "/picture"))
    }
    if (!dir.exists(paste0(task_path, "/output", "/store"))) {
        dir.create(paste0(task_path, "/output", "/store"))
    }
    if (!dir.exists(paste0(task_path, "/output", "/store", "/", prefix, "_", type))) {
        dir.create(paste0(task_path, "/output", "/store", "/", prefix, "_", type))
    }

    if (generate == F) {
        gene2cnv <- rownames(object) %>%
            AnnoProbe::annoGene(., "SYMBOL", species) %>%
            arrange(chr, start) %>%
            select(SYMBOL, chr, start, end) %>%
            distinct(SYMBOL, .keep_all = T)
        write.table(gene2cnv,
            file = paste0(task_path, "/output", "/store", "/", prefix, "_gene_file.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F
        )
        exp2cnv <- FetchData(object,
            vars = gene2cnv$SYMBOL,
            cells = colnames(object),
            layer = "count"
        ) %>% t()
        write.table(exp2cnv,
            file = paste0(task_path, "/output", "/store", "/", prefix, "_exp_file.txt"),
            sep = "\t",
        )
        anno2cnv <- object@meta.data %>%
            mutate(cellnames2cnv = rownames(.)) %>%
            select(cellnames2cnv, !!sym(celltype_col))
        write.table(anno2cnv,
            file = paste0(task_path, "/output", "/store", "/", prefix, "_anno_file.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F
        )
        rm(list = c("exp2cnv", "gene2cnv", "anno2cnv"))
    }

    message("Succeed to prepare the input")
    # 2 run the program
    future::plan("multisession", workers = ncores)
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = paste0(task_path, "/output", "/store", "/", prefix, "_exp_file.txt"), # raw counts matrix
        annotations_file = paste0(task_path, "/output", "/store", "/", prefix, "_anno_file.txt"), # an annotations file which indicates which cells are tumor vs. normal.
        delim = "\t",
        gene_order_file = paste0(task_path, "/output", "/store", "/", prefix, "_gene_file.txt"), # gene/chromosome positions file
        ref_group_names = normal_cell # reference cell name
    )
    message("Succeed to Create Infercnv Object")
    if (type == "samples") {
        infercnv_obj <- infercnv::run(infercnv_obj,
            cutoff = 0.1,
            # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
            out_dir = paste0(task_path, "/output", "/store", "/", prefix, "_", type),
            cluster_by_groups = TRUE,
            analysis_mode = type,
            denoise = TRUE,
            HMM = TRUE
        )
    } else {
        infercnv_obj <- infercnv::run(infercnv_obj,
            cutoff = 0.1,
            out_dir = paste0(task_path, "/output", "/store", "/", prefix, "_", type),
            cluster_by_groups = F,
            analysis_mode = type,
            hclust_method = "ward.D2",
            tumor_subcluster_partition_method = "random_trees",
            tumor_subcluster_pval = 0.05,
            denoise = TRUE,
            HMM = TRUE
        )
    }
    saveRDS(infercnv_obj, file = paste0(task_path, "/output", "/store", "/", prefix, "_infercnv_obj.RDS"))
    message("Succeed to run Infercnv")
    infercnv::plot_cnv(infercnv_obj,
        out_dir = paste0(task_path, "/output", "/picture", "/", prefix, "_", type),
        output_filename = paste0(prefix, "_infercnv.plot"),
        plot_chr_scale = F, # 是否画染色体全长
        output_format = "pdf",
        x.range = "auto",
        # x.center = 1,
        title = "infercnv",
        color_safe_pal = FALSE
    )
    return(infercnv_obj)
}

#' sc_cnv_copykat is a function to call copykat package to call Copy-number variants from
#' single cell data. It requires a Seurat object, an output path,
#' cell type column name, normal cell type name for normalization and number
#' of cores to use. It does preprocessing, then runs the copykat pipeline and generate
#' basic plots to visualize the results.
#'
#' @param object an object of class Seurat
#' @param task_path the output path
#' @param celltype_col the cell type column name
#' @param normal_cell the normal cell type names
#' @param ncores number of cores
#'
#' @return copykat.test object
sc_cnv_copykat <- function(object,
                           task_path,
                           prefix,
                           celltype_col,
                           normal_cell,
                           ncores) {
    # 0 check the setting
    library(copykat)
    assertthat::assert_that(class(object) == "Seurat")
    species <- check_species(object)
    assertthat::assert_that(species %in% c("human", "mouse"))
    assertthat::assert_that(all(celltype_col %in% colnames(object@meta.data)))
    assertthat::assert_that(all(normal_cell %in% unique(object@meta.data[, celltype_col])))
    ncores <- ncores %||% 20
    # 1 prepare the input
    if (!dir.exists(task_path)) {
        dir.create(task_path)
    }
    if (!dir.exists(paste0(task_path, "/input"))) {
        dir.create(paste0(task_path, "/input"))
    }
    if (!dir.exists(paste0(task_path, "/output"))) {
        dir.create(paste0(task_path, "/output"))
    }
    if (!dir.exists(paste0(task_path, "/output", "/picture"))) {
        dir.create(paste0(task_path, "/output", "/picture"))
    }
    if (!dir.exists(paste0(task_path, "/output", "/store"))) {
        dir.create(paste0(task_path, "/output", "/store"))
    }
    gene2cnv <- rownames(object) %>%
        AnnoProbe::annoGene(., "SYMBOL", species) %>%
        arrange(chr, start) %>%
        select(SYMBOL, chr, start, end) %>%
        distinct(SYMBOL, .keep_all = T)
    exp2cnv <- FetchData(object,
        vars = gene2cnv$SYMBOL,
        cells = colnames(object),
        layer = "count"
    ) %>% t()
    write.table(exp2cnv,
        file = paste0(task_path, "/output", "/store", "/", prefix, "_exp_file.txt"),
        sep = "\t",
        quote = F,
        row.names = T,
        col.names = T
    )
    # exp2cnv <- read.table(file = paste0(task_path, "/output", "/store","/", prefix, "_exp_file.txt"), sep = "\t", row.names = T, col.names = T)
    anno2cnv <- object@meta.data %>%
        mutate(cellnames2cnv = rownames(.)) %>%
        dplyr::select(cellnames2cnv, !!sym(celltype_col))
    normal2cnv <- anno2cnv %>%
        dplyr::filter(!!sym(celltype_col) %in% normal_cell) %>%
        dplyr::pull(cellnames2cnv)
    assertthat::assert_that(length(normal2cnv) > 50)
    message("Succeed to prepare for copykat")
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
    saveRDS(copykat_obj, file = paste0(task_path, "/output", "/store", "/", prefix, "_copykat_obj.RDS"))
    message("Succeed to run copykat")
    # 肿瘤是否恶性
    pred.test <- data.frame(copykat.test$prediction)
    pred.test <- lefj_join(
        x = pred.test, y = anno2cnv,
        dplyr::join_by(cell.names == cellnames2cnv)
    ) # copykat会自动过滤基因数少于100的细胞
    saveRDS(pred.test, file = paste0(task_path, "/output", "/store", "/", prefix, "_pred.test.RDS"))
    print(table(pred.test[, c(celltype_col, "copykat.pred")]))
    # cnv事件坐标与检测量
    CNA.test <- data.frame(copykat.test$CNAmat)
    saveRDS(CNA.test, file = paste0(task_path, "/output", "/store", "/", prefix, "_CNA.test.RDS"))
    message("Succeed to analysis of groups for copykat")
    # 3 plot the result
    ## chromosome color
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
    chr <- as.numeric(CNA.test$chrom) %% 2 + 1
    rbPal1 <- colorRampPalette(c("black", "grey"))
    CHR <- rbPal1(2)[as.numeric(chr)]
    chr1 <- cbind(CHR, CHR)
    ## scale color
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    com.preN <- pred.test$copykat.pred
    pred <- rbPal5(2)[as.numeric(factor(com.preN))]
    cells <- rbind(pred, pred)
    col_breaks <- c(
        seq(-1, -0.4, length = 50),
        seq(-0.4, -0.2, length = 150),
        seq(-0.2, 0.2, length = 600),
        seq(0.2, 0.4, length = 150),
        seq(0.4, 1, length = 50)
    )
    pdf(paste0(task_path, "/output", "/picture", "/", prefix, "_CNA_tumor_normal.pdf"),
        width = 14, height = 16
    )
    heatmap.3(t(CNA.test[, 4:ncol(CNA.test)]),
        dendrogram = "r",
        distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"),
        hclustfun = function(x) hclust(x, method = "ward.D2"),
        ColSideColors = chr1,
        RowSideColors = cells,
        Colv = NA, Rowv = TRUE,
        notecol = "black",
        col = my_palette,
        breaks = col_breaks,
        key = TRUE,
        keysize = 1, density.info = "none", trace = "none",
        cexRow = 0.1, cexCol = 0.1, cex.main = 1, cex.lab = 0.1,
        symm = F, symkey = F, symbreaks = T, cex = 1, cex.main = 4, margins = c(10, 10)
    )
    legend("topright",
        paste("pred.", names(table(com.preN)), sep = ""),
        pch = 15,
        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1],
        cex = 0.6, bty = "n"
    )
    dev.off()
    # 肿瘤根据层次聚类分亚群
    tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred == "aneuploid")]
    tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
    hcc <- hclust(parallelDist::parDist(t(tumor.mat), threads = 4, method = "euclidean"),
        method = "ward.D2"
    )
    hc.umap <- cutree(hcc, 2)
    saveRDS(list(hclust = hcc, hcut = hc.umap), file = paste0(task_path, "/output", "/store", "/", prefix, "_CNA_hclust_result.RDS"))
    message("Succeed to analysis of subclones for copykat")
    rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
    subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
    cells <- rbind(subpop, subpop)
    pdf(paste0(task_path, "/output", "/picture", "/", prefix, "_CNA_tumor_subclone.pdf"),
        width = 14, height = 16
    )
    heatmap.3(t(tumor.mat),
        dendrogram = "r",
        distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"),
        hclustfun = function(x) hclust(x, method = "ward.D2"),
        ColSideColors = chr1, RowSideColors = cells, Colv = NA, Rowv = TRUE,
        notecol = "black", col = my_palette, breaks = col_breaks, key = TRUE,
        keysize = 1, density.info = "none", trace = "none",
        cexRow = 0.1, cexCol = 0.1, cex.main = 1, cex.lab = 0.1,
        symm = F, symkey = F, symbreaks = T, cex = 1, cex.main = 4, margins = c(10, 10)
    )
    legend("topright", c("c1", "c2"),
        pch = 15,
        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4],
        cex = 0.9, bty = "n"
    )
    dev.off()
    return(copykat.test)
}




sc_cnv_score <- function(expr, methods) {
    # https://zhuanlan.zhihu.com/p/433234064
    # https://www.jianshu.com/p/1fa1fd4f97ff
    expr %<>% as.matrix()
    if (methods == 1) {
        expr.scale <- scale(t(expr))
        tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min), "-")
        tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale, 2, min)
        expr_1 <- t(2 * sweep(tmp1, 2, tmp2, "/") - 1)
        cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
        colnames(cnv_score) <- "cnv_score"
        cnv_score <- rownames_to_column(cnv_score, var = "cell")
        return(cnv_score)
    } else if (methods == 2) {
        expr[expr > 0 & expr < 0.3] <- 2 # complete loss. 2pts
        expr[expr >= 0.3 & expr < 0.7] <- 1 # loss of one copy. 1pts
        expr[expr >= 0.7 & expr < 1.3] <- 0 # Neutral. 0pts
        expr[expr >= 1.3 & expr <= 1.5] <- 1 # addition of one copy. 1pts
        expr[expr > 1.5 & expr <= 2] <- 2 # addition of two copies. 2pts
        expr[expr > 2] <- 2 # addition of more than two copies. 2pts
        cnv_score <- as.data.frame(colSums(expr))
        colnames(cnv_score) <- "cnv_score"
        cnv_score$cell <- rownames(cnv_score)
        cnv_score <- cnv_score[, c("cell", "cnv_score")]
        return(cnv_score)
    }
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
