sc_cnv <- function(object, method, output) {
    assertthat::assert_that(class(object) == "Seurat")

    sc_cnv_copykat <- function() {
        library(copykat)
    }
}



sc_cnv_infercnv <- function(outpath,
                            celltype_col,
                            species,
                            exp_file,
                            annotation_file,
                            gene_file,
                            denoise = T,
                            HMM = T) {
    assertthat::assert_that(species %in% c("human", "mouse"))
    assertthat::assert_that(all(celltype_col %in% colnames(object@meta.data)))
    dir.create(paste0(outpath, "/input"))
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
    dir.create(paste0(outpath, "/output"))
    dir.create(paste0(outpath, "/output", "/picture"))
    dir.create(paste0(outpath, "/output", "/samples"))
    dir.create(paste0(outpath, "/output", "/subclusters"))
    library(infercnv)
    library(ggplot2)
    future::plan("multiprocess", workers = 20)
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = exp_file, # raw counts matrix
        annotations_file = annotation_file, # an annotations file which indicates which cells are tumor vs. normal.
        delim = "\t",
        gene_order_file = gene_file, # gene/chromosome positions file
        ref_group_names = c("Microglia/Macrophage", "Oligodendrocytes (non-malignant)") # 作为参考的细胞名
    )

    infercnv_obj <- infercnv::run(infercnv_obj,
        cutoff = 0.1,
        # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
        out_dir = paste0(outpath, "/output", "/samples"),
        cluster_by_groups = TRUE,
        analysis_mode = "samples",
        denoise = TRUE,
        HMM = TRUE
    )
    infercnv_obj <- infercnv::run(infercnv_obj,
        cutoff = 0.1,
        out_dir = paste0(outpath, "/output", "/subclusters"),
        cluster_by_groups = F,
        analysis_mode = "subclusters",
        hclust_method = "ward.D2",
        tumor_subcluster_partition_method = "random_trees",
        tumor_subcluster_pval = 0.05,
        denoise = TRUE,
        HMM = TRUE
    )
    save(infercnv_obj, file = paste0(outpath, "/output", "/picture", "/infercnv_obj.RData"))
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

# https://zhuanlan.zhihu.com/p/433234064
# https://www.jianshu.com/p/1fa1fd4f97ff
sc_cnv_score <- function() {

}
