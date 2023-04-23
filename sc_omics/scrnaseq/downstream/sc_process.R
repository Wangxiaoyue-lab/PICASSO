#' Process and read data from file
#'
#' @param filename The name of the file to read
#' @param type The type of file to read
#' @param project The project name for the Seurat object
#' @param min.cells Minimum number of cells for creating Seurat object
#'
#' @return A Seurat object
process_read <- function(filename,
                         type,
                         project,
                         min.cells = 3, ...) {
    read_h5 <- function(filename) {
        Seurat::Read10X_h5(filename = filename, use.names = T) %>%
            CreateSeuratObject(project = project, min.cells = min.cells, ...)
    }
    read_h5ad <- function(filename) {
        library(SeuratDisk)
        Convert(filename, "h5seurat", overwrite = T, assay = "RNA")
        LoadH5Seurat(filename %>%
            str_sub(., 1, nchar(.) - 2) %>%
            str_c(., "seurat"))
    }
    read_10x <- function(filename) {
        if (length(list.files(filename)) < 3) {
            stop("The 10x files is not right")
        }
        Read10X(filename) %>%
            CreateSeuratObject(project = project, min.cells = min.cells, ...)
    }
    read_loom <- function(filename) {
        library(SeuratDisk)
        library(SeuratObject)
        Connect(filename, mode = "r+") %>%
            .[[matrix]] %>%
            as.Seurat()
    }
    # read_mtx <- function(filename){
    #    readMM(filename)
    # }
    read_table <- function(filename) {
        library(Matrix)
        as.matrix(filename) %>%
            as(., "dgCMatrix") %>%
            CreateSeuratObject(project = project, min.cells = min.cells, ...)
    }
    read_seurat <- switch(type,
        "h5" = read_h5,
        "h5ad" = read_h5ad,
        "10x" = read_10x,
        "loom" = read_loom,
        #' mtx'=read_mtx,
        "table" = read_table
    )
    read_seurat(filename)
}


process_read_v5 <- function(filename,
                            type,
                            project,
                            min.cells = 3,
                            store_path, ...) {
    library(BPCells)
    library(Seurat)
    library(SeuratObject)
    library(SeuratDisk)
    options(Seurat.object.assay.version = "v5")
    read_h5 <- function(filename){
        mat <- open_matrix_10x_hdf5(path = filename) %>% 
            write_matrix_dir(mat = ., dir = store_path, overwrite = TRUE)
    }
    read_10x <- function(filename){
        mat <- open_matrix_10x(path = filename) %>% 
            write_matrix_dir(mat = ., dir = store_path, overwrite = TRUE)
    }
    object <- open_matrix_dir(dir = store_path) %>% 
        CreateSeuratObject(counts = . )
}

process_to_v5 <- function(object,store_path) {
    options(Seurat.object.assay.version = "v5")
    if(!dir.exists(store_path)){
        dir.create(store_path, recursive = T)
    }
    write_matrix_dir(mat = object[["RNA"]]$counts, 
        dir = store_path)
    meta.data <- object@meta.data
    counts.mat <- open_matrix_dir(dir = store_path)
    object <- CreateSeuratObject(counts = counts.mat)
    object <- AddMetaData(object,meta.data)
    gc()
    return(object)
}


#' Process Seurat object to 3 files
#'
#' @param object The Seurat object to process
#' @param output The output directory for the 3 files
#'
#' @return 3 files for the Seurat object
process_to_3files <- function(object,
                              output) {
    if (!requireNamespace("DropletUtils", quietly = TRUE)) {
        BiocManager::install("DropletUtils")
    }
    library(DropletUtils)
    write10xCounts(output, object[["RNA"]]@counts, version = 3)
}


#' Process Seurat object
#'
#' @param object The Seurat object to process
#' @param npcs Number of principal components to use
#' @param resolutions Resolutions for clustering
#' @param future Logical, whether to use future package for parallel processing
#' @param run_harmony Logical, whether to run Harmony for batch correction
#' @param run_sctransform Logical, whether to run SCTransform for normalization
#' @param group_in_harmony Grouping variable for Harmony batch correction
#' @param vars_to_regress Variables to regress out during normalization
#' @param verbose Logical, whether to print progress messages
#'
#' @return A processed Seurat object
process_process <- function(object,
                            npcs = 20,
                            resolutions = c(0.1, 0.2, 0.3, 0.5),
                            future = FALSE,
                            run_harmony = TRUE,
                            run_sctransform = TRUE,
                            group_in_harmony = "orig.ident",
                            vars_to_regress = c("percent_mito", "S.Score", "G2M.Score"),
                            verbose = F, ...) {
    if (future) {
        options(future.globals.maxSize = check_size_future(object))
        plan("multisession")
    }
    assay_use <<- switch(run_sctransform + 1,
        "RNA",
        "SCT"
    )
    ident_value <- paste0(assay_use, "_snn_res.", min(resolutions))
    if (run_sctransform == TRUE) {
        object <- SCTransform(object, vars.to.regress = vars_to_regress, verbose = verbose) %>%
            RunPCA(verbose = verbose)
    } else {
        object <- NormalizeData(object, verbose = verbose) %>%
            FindVariableFeatures(verbose = verbose) %>%
            ScaleData(features = row.names(object), vars.to.regress = vars_to_regress, verbose = verbose) %>%
            RunPCA(verbose = verbose)
    }

    reduction_use <- switch(run_harmony + 1,
        "pca",
        "harmony"
    )

    if (run_harmony) {
        library(harmony)
        object <- object %>%
            RunHarmony(group.by.vars = group_in_harmony, dims.use = 1:npcs, assay.use = assay_use)
    }
    object <- RunUMAP(object, reduction = reduction_use, dims = 1:npcs, verbose = verbose) %>%
        FindNeighbors(reduction = reduction_use, dims = 1:npcs, verbose = verbose) %>%
        FindClusters(resolution = resolutions, verbose = verbose) %>%
        SetIdent(value = ident_value)

    return(object)
}



#' Add metadata to Seurat object
#'
#' @param object The Seurat object to add metadata to
#' @param new.meta New metadata to add
#' @param by.o Old column name for joining
#' @param by.n New column name for joining
#' @param type Type of metadata to add, either "sample" or "cell"
#' @param filter Logical, whether to filter cells not present in new metadata
#'
#' @return A Seurat object with added metadata
process_add_meta.data <- function(object,
                                  new.meta,
                                  by.o = NULL, # old/object
                                  by.n, # new
                                  type = c("sample", "cell"),
                                  filter = FALSE) {
    join_ <- ifelse(filter, inner_join, left_join)
    object@meta.data$cell_names <- row.names(object@meta.data)
    if (type == "cell") {
        by.o <- "cell_names"
    }
    meta.filt <- object@meta.data %>%
        join_(., new.meta, join_by(!!sym(by.o) == !!sym(by.n)))
    row.names(meta.filt) <- meta.filt$cell_names
    object <- subset(object, subset = cell_names %in% meta.filt$cell_names)
    object[[colnames(meta.filt)]] <- meta.filt
    # object <- AddMetaData(object, metadata = meta.filt)
    return(object)
}

#' Annotate Seurat object
#'
#' @param object The Seurat object to annotate
#' @param col.id Column name for old identities
#' @param col.new Column name for new identities
#' @param split Logical, whether to split object by new identities
#'
#' @return An annotated Seurat object or a list of split Seurat objects
process_annotation <- function(object,
                               col.id,
                               col.new,
                               split = F, ...) {
    ident.pairs <- tryCatch(
        expr = as.list(x = ...),
        error = function(e) {
            return(list(...))
        }
    ) %>% unlist()
    col.old <- object@meta.data %>%
        pull(col.id) %>%
        unique()
    new_minus_old <- setdiff(names(ident.pairs), col.old)
    if (length(new_minus_old) > 0) {
        {
            stop(paste(
                "The identities",
                paste0(new_minus_old, collapse = ","), "do not exist"
            ))
            return(list(...))
        } %>% unlist()
    }
    col.old <- object@meta.data %>%
        pull(!!sym(col.id)) %>%
        unique()
    new_minus_old <- setdiff(names(ident.pairs), col.old)
    if (length(new_minus_old) > 0) {
        stop(paste(
            "The identities",
            paste0(new_minus_old, collapse = ","), "do not exist"
        ))
    }
    old_minus_new <- setdiff(col.old, names(ident.pairs))
    if (length(old_minus_new) > 0) {
        ident.pairs[[old_minus_new]] <- "unknown"
        warning(paste(
            "The identities",
            paste0(old_minus_new, collapse = ","), "do not provide the clear celltype"
        ))
    }
    object@meta.data %<>% mutate(
        !!col.new := map(!!sym(col.id), ~ {
            ident.pairs[as.character(.x)]
        })
    )

    if (split == F) {
        return(object)
    } else {
        return(SplitObject(object, split.by = col.new))
    }
}


# process_integration

#' Find marker genes in Seurat object
#'
#' @param object The Seurat object to find marker genes in
#' @param is_all Logical, whether to find all markers or markers for specific identities
#' @param future Logical, whether to use future package for parallel processing
#' @param ident Identity to use for finding all markers
#' @param loop_var Loop variable for comparing conditions within identities
#' @param species Species for adding gene annotations
#'
#' @return A data frame of marker genes
process_find_markers <- function(object,
                                 is_all = TRUE, # whether FindMarkers or FindAllMarkers
                                 future = FALSE, # whether use future
                                 ident = NULL, # the ident for FindAllMarkers
                                 loop_var = NULL, # each celltype or other ident
                                 species, # the species
                                 ...) {
    # @is_all：FindAllMarkers or FindMarkers
    # @loop_var: compare conditions within some idents

    if (future) {
        options(future.globals.maxSize = check_size_future(object))
        plan("multicore")
    }
    assay_use <<- ifelse(exists("assay_use"), assay_use %||%
        DefaultAssay(object), DefaultAssay(object))

    add_fun <- function(marker_genes) {
        annotations <- utils_read_refdata(species, "annotations")
        # Add a column with the percentage difference
        marker_genes <- Add_Pct_Diff(marker_genes) %>%
            # Create a new column with the type of marker gene
            mutate(type = if_else(avg_log2FC > 0, "up", "down")) %>%
            # Use the left_join() function to add the gene name and description
            left_join(
                unique(annotations[, c("gene_name", "description")]),
                by = c("gene" = "gene_name")
            )
        return(marker_genes)
    }
    if (assay_use == "SCT") {
        object <- PrepSCTFindMarkers(object)
    }
    if (!is.null(ident)) {
        object <- SetIdent(object, value = ident)
    }
    if (is_all == TRUE) {
        # If all is TRUE, find all markers
        marker_genes <- FindAllMarkers(object, assay = assay_use, ...) %>% add_fun()
    } else {
        if (!is.null(loop_var)) {
            # If all is FALSE and loop_var is not NULL, find markers for each cluster
            data_list <- lapply(loop_var, function(x) {
                marker_genes <- FindMarkers(object, assay = assay_use, subset.ident = x, ...) %>%
                    rownames_to_column("gene") %>%
                    mutate(cluster = rep(x, nrow(.))) %>%
                    add_fun()
            })
            marker_genes <- do.call(rbind, data_list)
        } else {
            # If all is FALSE and loop_var is NULL, find markers
            marker_genes <- FindMarkers(object, assay = assay_use, ...) %>%
                rownames_to_column("gene") %>%
                add_fun()
        }
    }
    return(marker_genes)
}


# process_celltype



#' Process individual degrees of freedom
#'
#' @param object An object containing data to be processed
#' @param methods A character vector specifying the methods to use for processing
#' @param sample_id A character vector specifying the sample IDs
#' @param group_id A character vector specifying the group IDs
#' @param covariables A character vector specifying the covariables to consider (default: NULL)  For: ideas
#' @param need_scale A logical value indicating whether continuous variables need to be scaled (default: NULL) For: ideas
#' @param cell_select A character vector specifying the cells to select (default: NULL) For: deseq2
#' @param cluster_id A character vector specifying the cluster IDs (default: NULL) For: deseq2
#' @param filepath A character vector specifying the file path (default: NULL) For: deseq2
#'
#' @return A list containing the results of processing individual degrees of freedom
#'
#' @export
process_individual_deg <- function(object,
                                   methods,
                                   sample_id,
                                   group_id,
                                   covariables = NULL, # ideas
                                   need_scale = NULL, # ideas
                                   cell_select = NULL, # deseq2
                                   cluster_id = NULL, # deseq2
                                   filepath = NULL # deseq2
) {
    process_ideas <- function(object,
                              sample_id, # the id of patients of samples
                              group_id, # the group of patients(i.e. tumor or drug)
                              covariables, # the covariables should be considered(i.e. gender or age)
                              need_scale # the continous variables need to be scaled
    ) {
        library(ideas)
        library(SingleCellExperiment)
        library(doParallel)
        library(doRNG)
        library(foreach)
        ncores <- 20
        registerDoParallel(cores = ncores)
        options(mc.core = ncores)
        RNGkind("L'Ecuyer-CMRG")
        var2test <- "diagnosis" # 检验变量
        var2adjust <- covariables # 协变量
        var2test_type <- "binary" # 检验类型
        var_per_cell <- c("rd") # 纠正深度
        # 表达矩阵
        count_matrix <- object[["RNA"]]@counts %>%
            as.matrix() %>%
            as.data.frame() %>%
            rowwise() %>%
            mutate(sum = sum(c_across(everything()) == 0)) %>%
            filter(sum < 0.9 * (ncol(.) - 1)) %>%
            select(-sum) %>%
            as.data.frame() %>%
            set_rownames(row.names(object[["RNA"]]@counts)) %>%
            as.matrix()
        # 细胞元数据
        meta_cell <- object@meta.data %>%
            mutate(
                cell_id = rownames(.),
                individual = !!sym(sample_id),
                diagnosis = !!sym(group_id),
                rd = colSums(count_matrix)
            )
        # 病人元数据
        meta_ind <- meta_cell %>%
            select(individual, diagnosis, any_of(covariables)) %>%
            distinct() %>%
            mutate(across(need_scale, scale))
        dist1 <- ideas_dist(count_matrix,
            meta_cell,
            meta_ind,
            var_per_cell,
            var2test,
            var2test_type,
            d_metric = "Was",
            fit_method = "nb"
        )
        pval_ideas <- permanova(dist1,
            meta_ind,
            var2test,
            var2adjust,
            var2test_type,
            n_perm = 999,
            r.seed = 2023
        )
        p_ideas_res <- data.frame(
            genes = row.names(count_matrix),
            pval = pval_ideas
        ) %>%
            filter(!is.na(pval)) %>%
            arrange(desc(pval)) %>%
            filter(pval < 0.05)
        return(p_ideas_res)
    }
    process_deseq2 <- function(object,
                               cell_select,
                               cluster_id, # the cluster repesented celltype
                               sample_id, # the id of patients of samples
                               group_id, # the group of patients like tumor or drug
                               covariables, # the covariables should be considered like gender or age
                               need_scale, # the continous variables need to be scaled
                               filepath) {
        library(Matrix.utils)
        library(Matrix)
        library(SingleCellExperiment)
        library(apeglm)
        library(png)
        sce <- SingleCellExperiment(
            assays = list(counts = object@assays$RNA@counts),
            colData = object@meta.data %>%
                mutate(
                    sample_id = !!sym(sample_id),
                    group_id = !!sym(group_id),
                    cluster_id = !!sym(cluster_id)
                )
        )
        groups <- colData(sce)[, c("cluster_id", "sample_id")]

        # 生成样本级的元数据
        sce$cluster_id <- factor(sce$cluster_id)
        kids <- purrr::set_names(levels(sce$cluster_id))
        nk <- length(kids)
        sce$sample_id <- factor(sce$sample_id)
        sids <- purrr::set_names(levels(sce$sample_id))
        ns <- length(sids)
        m <- match(sids, sce$sample_id)
        ei <- data.frame(colData(sce)[m, ],
            as.numeric(table(sce$sample_id)),
            row.names = NULL
        ) %>%
            select(-"cluster_id")

        # Aggregate across cluster-sample groups
        pb <- aggregate.Matrix(t(counts(sce)),
            groupings = colData(sce)[, c("cluster_id", "sample_id")],
            fun = "mean"
        )
        splitf <- sapply(
            stringr::str_split(rownames(pb),
                pattern = "_",
                n = 2
            ),
            `[`, 1
        )
        pb <- split.data.frame(
            pb,
            factor(splitf)
        ) %>%
            lapply(function(u) {
                set_colnames(
                    t(u),
                    stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")
                )
            })

        # deg 元数据
        get_sample_ids <- function(x) {
            pb[[x]] %>%
                colnames()
        }
        de_samples <- map(1:length(kids), get_sample_ids) %>%
            unlist()
        samples_list <- map(1:length(kids), get_sample_ids)
        get_cluster_ids <- function(x) {
            rep(names(pb)[x],
                each = length(samples_list[[x]])
            )
        }
        de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
            unlist()
        gg_df <- data.frame(
            cluster_id = de_cluster_ids,
            sample_id = de_samples
        )
        gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")])
        metadata <- gg_df %>%
            dplyr::select(cluster_id, sample_id, group_id)

        ## 提取某一亚群细胞
        clusters <- levels(metadata$cluster_id)
        cluster_x <- cell_select %||% clusters[1]
        cluster_metadata <- metadata[which(metadata$cluster_id == cluster_x), ]
        rownames(cluster_metadata) <- cluster_metadata$sample_id
        counts <- pb[[cluster_x]]
        cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
        all(rownames(cluster_metadata) == colnames(cluster_counts))
        dds <- DESeqDataSetFromMatrix(cluster_counts,
            colData = cluster_metadata,
            design = ~group_id
        )
        rld <- rlog(dds, blind = TRUE)
        pdf(paste0(filepath))
        DESeq2::plotPCA(rld, intgroup = "group_id")
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        pheatmap(rld_cor,
            annotation = cluster_metadata[, c("group_id"), drop = F]
        )
        dev.off()
        dds <- DESeq(dds)
        cluster_metadata$group_id <- factor(cluster_metadata$group_id)
        contrast <- c(
            "group_id",
            levels(cluster_metadata$group_id)[2],
            levels(cluster_metadata$group_id)[1]
        )
        res <- results(dds,
            contrast = contrast,
            alpha = 0.05
        )
        res <- lfcShrink(dds,
            contrast = contrast,
            type = "ashr"
        )
        return(res)
    }
    process_mast <- function(...) {
        stop("the method are planning to be design")
    }
    deg_function <- switch(methods,
        "ideas" = process_ideas,
        "pseudobulk" = process_deseq2,
        "mast" = process_mast
    )
}
