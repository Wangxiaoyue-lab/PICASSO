sc_score_function <- function(...) {
    cat("# function list\n")
    cat("#-- sc_score_seurat: scoring by seurat\n")
    cat("#-- sc_score_aucell: scoring by aucell\n")
    cat("#-- sc_score_ucell: scoring by ucell\n")
    cat("#-- sc_score_pagoda2: scoring by pagoda2\n")
    cat("#-- sc_score_maya: scoring by maya\n")
    cat("#-- sc_score_jasmine: scoring by jasmine\n")
    cat("#-- warning: The following methods are slow!")
}

# Seurat - AddModuleScore
sc_score_seurat <- function(object, genes_list) {
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 2 run
    object <- Seurat::AddModuleScore(object,
        features = genes_list,
        ctrl = 100,
        name = paste0(names(genes_list), "_seurat")
    )
    return(object)
}

# AUCell
sc_score_aucell <- function(object, genes_list) {
    library(AUCell)
    set.seed(1)
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 2 run
    matrix <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "data"
    ) %>% t()
    cells_rankings <- AUCell_buildRankings(matrix, nCores = 1, plotStats = TRUE)
    lapply(seq_along(genes_list), function(g_l) {
        cells_AUC <- AUCell_calcAUC(genes_list[[g_l]], cells_rankings)
        signature_exp <- data.frame(t(getAUC(cells_AUC)))
        object@meta.data[[paste0(names(genes_list)[g_l], "_aucell")]] <- signature_exp[, 1]
    })
}

# UCell
sc_score_ucell <- function(object, genes_list) {
    library(UCell)
    library(GSEABase)
    library(BiocParallel)
    set.seed(1)
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 2 run
    object <- AddModuleScore_UCell(object,
        features = genes_list,
        assay = NULL,
        slot = "data",
        BPPARAM = NULL,
        ncores = 1,
        chunk.size = 1000
    )
    return(object)
}

# singscore
sc_score_singscore <- function(object, genes_list) {
    library(singscore)
    set.seed(1)
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    matrix <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "data"
    ) %>% t()
    singscore.rank <- singscore::rankGenes(as.data.frame(matrix))
    singscore.Score <- lapply(seq_along(genes_list), function(g_l) {
        singscore.scores_g_l <- singscore::simpleScore(singscore.rank,
            upSet = genes_list[[g_l]],
            centerScore = F
        ) %>%
            as.data.frame() %>%
            dplyr::select(TotalScore)
    }) %>% do.call(cbind, .)
    colnames(singscore.Score) <- paste0(names(genes_list), "_singscore")
    object@meta.data <- cbind(object@meta.data, singscore.Score)
    return(object)
}


sc_score_plage <- function(...) {
    # too slow
    # 0 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    n_cores <- n_cores %||% 4
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 1 extract data
    counts <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()
    # 2 run
    score <- GSVA::gsva(counts,
        genes_list,
        method = "plage",
        parallel.sz = n_cores,
        verbose = T
    ) %>% as.data.frame()
    colnames(score) <- paste0(colnames(score), "_gsva")
    object@meta.data <- cbind(object@meta.data, score)
    return(object)
}
sc_score_zscore <- function(...) {
    # too slow
    # 0 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    n_cores <- n_cores %||% 4
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 1 extract data
    counts <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()
    # 2 run
    score <- GSVA::gsva(counts,
        genes_list,
        method = "zscore",
        parallel.sz = n_cores,
        verbose = T
    ) %>% as.data.frame()
    colnames(score) <- paste0(colnames(score), "_gsva")
    object@meta.data <- cbind(object@meta.data, score)
    return(object)
}
sc_score_gsea <- function(...) {
    # too slow
}
sc_score_ssgsea <- function(...) {
    # too slow
    # 0 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    n_cores <- n_cores %||% 4
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 1 extract data
    counts <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()
    # 2 run
    score <- GSVA::gsva(counts,
        genes_list,
        method = "ssgsea",
        parallel.sz = n_cores,
        verbose = T
    ) %>% as.data.frame()
    colnames(score) <- paste0(colnames(score), "_gsva")
    object@meta.data <- cbind(object@meta.data, score)
    return(object)
}
sc_score_gsva <- function(object, genes_list, n_cores = NULL) {
    # too slow
    # 0 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    n_cores <- n_cores %||% 4
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 1 extract data
    counts <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()
    # 2 run
    score <- GSVA::gsva(counts,
        genes_list,
        method = "gsva",
        parallel.sz = n_cores,
        verbose = T
    ) %>% as.data.frame()
    colnames(score) <- paste0(colnames(score), "_gsva")
    object@meta.data <- cbind(object@meta.data, score)
    return(object)
}

# The function should be checked
# sc_score_vision <- function(...) {
# too slow
# 0 check
#    assertthat::assert_that(class(object) == "Seurat")
#    assertthat::assert_that(class(genes_list) == "list")
#    n_cores <- n_cores %||% 4
#    genes_list <- lapply(genes_list, function(g) {
#        intersect(genes_list[[g]], rownames(object))
#    }) %>% list_clean()
#    # 1 extract data
#    counts <- FetchData(object,
#        vars = rownames(object),
#        cells = colnames(object),
#        layer = "counts"
#    ) %>% t()
#    # 2 run
#    vis = VISION::Vision(counts,            ## Gene X Cell
#        # data.frame; sparseMatrix; dgeMatrix; ExpressionSet; SummarizedExperiment; Seurat
#        signatures = gSets_path,
#        projection_method = 'UMAP',
#        sig_gene_threshold=0)
#    options(mc.cores=n_cores)
#    vis = VISION::analyze(vis)
#    score = t(vis@SigScores)    ## pathway X cell
#    return(score)
# }


# Pagoda2
sc_score_pagoda2 <- function(object, genes_list, n_cores = NULL) {
    # https://github.com/sulab-wmu/PASBench/blob/master/R/tools.R
    library(pagoda2)
    # 0 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    n_cores <- n_cores %||% 4
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 1 extract data
    counts <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()
    # 2 construct the pagoda2 object
    nPcs <- min(round(ncol(counts) / 5), 5)
    # counts = apply(counts,2,function(x) {storage.mode(x) = 'integer'; x})
    p2 <- Pagoda2$new(counts, n.cores = n_cores, log.scale = F)
    p2$adjustVariance(plot = F)
    p2$calculatePcaReduction(nPcs = nPcs, use.odgenes = F, fastpath = F)
    # 3 run
    path_names <- c()
    env <- new.env(parent = globalenv())
    invisible(lapply(1:length(genes_list), function(i) {
        genes <- intersect(genes_list[[i]], rownames(counts))
        name <- paste0(names(genes_list[i]), i)
        if (length(genes) > 3) {
            assign(name, genes, envir = env)
            path_names <- c(path_names, name)
        }
    }))
    p2$testPathwayOverdispersion(
        setenv = env, verbose = T,
        recalculate.pca = T,
        min.pathway.size = 1
    )
    path_names <- names(p2@.xData$misc$pwpca)
    score <- matrix(NA, nrow = length(path_names), ncol = ncol(counts))
    rownames(score) <- path_names
    colnames(score) <- colnames(counts)
    for (i in 1:length(p2@.xData$misc$pwpca)) {
        if (!is.null(p2@.xData$misc$pwpca[[i]]$xp$score)) {
            score[i, ] <- as.numeric(p2@.xData$misc$pwpca[[i]]$xp$scores)
        }
    }
    object@meta.data <- cbind(object@meta.data, as.data.frame(score))
    return(object)
}

# MAYA
sc_score_maya <- function(object, genes_list = NULL, n_cores = NULL) {
    # https://github.com/One-Biosciences/MAYA/
    library(MAYA)
    # 0 modified the function
    study_pathways <- function(PCA_obj, compute_umap = T) {
        # compute activity mat
        activity_mat <- build_activity_mat(PCA_obj, scaled = F)

        # clustering and generate average score by matrix
        # knn.matrix <- RANN::nn2(t(activity_mat), t(activity_mat), k = 20, searchtype = "standard")[[1]]
        # jaccard.adj <- knn_jaccard(knn.matrix)
        # graph <- igraph::graph.adjacency(
        #    jaccard.adj,
        #    mode = "undirected",
        #    weighted = TRUE
        # )
        # set.seed(2016)
        # cluster_result <- leidenbase::leiden_find_partition(
        #    graph,
        #    verbose = FALSE,
        #    seed = 2016, num_iter = 2, partition_type = "ModularityVertexPartition", edge_weights = igraph::E(graph)$weight
        # )
        # cluster_result$membership <- paste0("C", cluster_result$membership)
        # tmp <- average_by_cluster(activity_mat, cluster_result$membership)

        # compute umap
        if (compute_umap) {
            umap <- run_umap(activity_mat)
        } else {
            umap <- NULL
        }
        return(list(cluster_matrix = NULL, clusters_annotation = NULL, activity_matrix = activity_mat, umap = umap))
        # return(list(cluster_matrix = tmp, clusters_annotation = cluster_result$membership, activity_matrix = activity_mat, umap = umap))
    }

    MAYA_pathway_analysis <- function(expr_mat, modules_list = NULL, min_cells_pct = 0.05, is_logcpm = T, nCores = 1, min_genes = 10, max_contrib = 0.5, compute_umap = T, scale_before_pca = T, all_PCs_in_range = F) {
        message("Running pathway analysis")

        #### Check parameters ####
        stopifnot(
            is.logical(is_logcpm),
            is.numeric(min_cells_pct), is.numeric(nCores)
        )

        # load modules list in Panglao if necessary
        if (is.null(modules_list)) {
            message("Loading HALLMARK from MSigDB")
            path <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "MAYA")
            modules_list <- read_gmt(path)
        }
        if (is.character(modules_list)) {
            if (modules_list == "hallmark") {
                message("Loading HALLMARK from MSigDB")
                path <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "MAYA")
                modules_list <- read_gmt(path)
            } else {
                if (modules_list == "kegg") {
                    message("Loading KEGG from MSigDB")
                    path <- system.file("extdata", "c2.cp.kegg.v7.4.symbols.gmt", package = "MAYA")
                    modules_list <- read_gmt(path)
                }
            }
        }
        # at this stage, it can only be a list
        stopifnot(is.list(modules_list))

        # run MAYA with pathways parameters
        suppressWarnings(PCA_obj <- run_activity_analysis(
            expr_mat = expr_mat,
            modules_list = modules_list,
            nb_comp_max = 5,
            min_cells_pct = min_cells_pct,
            min_module_size = min_genes,
            max_contrib = max_contrib,
            norm = !is_logcpm,
            nCores = nCores,
            scale_before_pca = scale_before_pca,
            all_PCs_in_range = all_PCs_in_range
        ))
        # analyze pathways
        if (length(PCA_obj) != 0) {
            annot <- study_pathways(PCA_obj, compute_umap = compute_umap)
        } else {
            annot <- NULL
        }

        return(c(annot, list(PCA_obj = PCA_obj)))
    }
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    n_cores <- n_cores %||% 4
    if (is.null(genes_list)) {
        genes_list <- "hallmark"
    } else {
        assertthat::assert_that(class(genes_list) == "list")
        genes_list <- lapply(genes_list, function(g) {
            intersect(genes_list[[g]], rownames(object))
        }) %>% list_clean()
    }
    # 2 extract data
    counts <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()
    # 3 run
    activity_summary <- MAYA_pathway_analysis(
        expr_mat = counts,
        modules_list = genes_list,
        is_logcpm = F,
        nCores = n_cores,
        compute_umap = F
    )
    score <- scale_0_1(activity_summary$activity_matrix) %>%
        t() %>%
        as.data.frame()
    colnames(score) <- paste0(colnames(score), "_maya")
    object@meta.data <- cbind(object@meta.data, score)
    return(object)
}

# JASMINE
sc_score_jasmine <- function(object, genes_list) {
    # 0 design function
    # https://github.com/NNoureen/JASMINE
    stringsAsFactors <- FALSE
    library(stringr)
    library(GSA)
    set.seed(1)
    ##################################### STRUCTURE of JASMINE ####################
    ### Function1:-  Calculating Mean Ranks for signature genes across each cell
    RankCalculation <- function(x, genes) {
        subdata <- x[x != 0] ### Removing Dropouts from single cell
        DataRanksUpdated <- rank(subdata) ### Calculating ranks of each signature gene per cell
        DataRanksSigGenes <- DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)] ### Shortling rank vector for signature genes
        CumSum <- ifelse(length(DataRanksSigGenes), mean(DataRanksSigGenes, na.rm = TRUE), 0) ### Calculating Mean of ranks for signature genes
        FinalRawRank <- CumSum / length(subdata) ### Normalizing Means by total coverage
        return(FinalRawRank)
    }
    #### Function2:- Calculating enrichment of signature genes across each cell 	(using odds ratio)
    ORCalculation <- function(data, genes) {
        GE <- data[which(rownames(data) %in% genes), ] ### Subsetting data for signature genes
        NGE <- data[-which(rownames(data) %in% genes), ] ### Subsetting data for non-signature genes
        SigGenesExp <- apply(GE, 2, function(x) length(x[x != 0])) ### Calculating Number of expressed Signature Genes per cell
        NSigGenesExp <- apply(NGE, 2, function(x) length(x[x != 0])) ### Calculating Number of expressed Non-Signature Genes per cell
        SigGenesNE <- nrow(GE) - SigGenesExp ### Calculating Number of Not expressed Signature Genes per cell
        SigGenesNE <- replace(SigGenesNE, SigGenesNE == 0, 1) ### Replacing Zero's with 1
        NSigGenesExp <- replace(NSigGenesExp, NSigGenesExp == 0, 1) ### Replacing Zero's with 1
        NSigGenesNE <- nrow(data) - (NSigGenesExp + SigGenesExp) ### Calculating Number of Not expressed Non-Signature Genes per cell
        NSigGenesNE <- NSigGenesNE - SigGenesNE
        OR <- (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp) ### Calculating Enrichment (Odds Ratio)
        return(OR)
    }
    #### Function3:- Calculating enrichment of signature genes across each cell (using Likelihood ratio)
    LikelihoodCalculation <- function(data, genes) {
        GE <- data[which(rownames(data) %in% genes), ]
        NGE <- data[-which(rownames(data) %in% genes), ]
        SigGenesExp <- apply(GE, 2, function(x) length(x[x != 0]))
        NSigGenesExp <- apply(NGE, 2, function(x) length(x[x != 0]))
        SigGenesNE <- nrow(GE) - SigGenesExp
        SigGenesNE <- replace(SigGenesNE, SigGenesNE == 0, 1)
        NSigGenesExp <- replace(NSigGenesExp, NSigGenesExp == 0, 1)
        NSigGenesNE <- nrow(data) - (NSigGenesExp + SigGenesExp)
        NSigGenesNE <- NSigGenesNE - SigGenesNE
        LR1 <- SigGenesExp * (NSigGenesExp + NSigGenesNE)
        LR2 <- NSigGenesExp * (SigGenesExp + SigGenesNE)
        LR <- LR1 / LR2
        return(LR)
    }
    ###  Function 4:- Scalar [0,1] Normalization of Means and Enrichment across set of cells
    NormalizationJAS <- function(JAS_Scores) {
        JAS_Scores <- (JAS_Scores - min(JAS_Scores)) / (max(JAS_Scores) - min(JAS_Scores))
        return(JAS_Scores)
    }
    ### Function 5:- Signature Scoring via JASMINE mergining Means and Enrichment
    JASMINE <- function(data, genes, method) {
        idx <- match(genes, rownames(data))
        idx <- idx[!is.na(idx)]
        if (length(idx) > 1) {
            RM <- apply(data, 2, function(x) RankCalculation(x, genes)) ### Mean RankCalculation for single cell data matrix
            RM <- NormalizationJAS(RM) ### Normalizing Mean Ranks

            if (method == "oddsratio") {
                OR <- ORCalculation(data, genes) ### Signature Enrichment Calculation for single cell data matrix (OR)
                OR <- NormalizationJAS(OR) ### Normalizing Enrichment Scores (OR)
                JAS_Scores <- (RM + OR) / 2
            } else if (method == "likelihood") {
                LR <- LikelihoodCalculation(data, genes) ### Signature Enrichment Calculation for single cell data matrix  (LR)
                LR <- NormalizationJAS(LR) ### Normalizing Enrichment Scores (LR)
                JAS_Scores <- (RM + LR) / 2
            }
            FinalScores <- data.frame(names(RM), JAS_Scores) ### JASMINE scores
            colnames(FinalScores)[1] <- "SampleID"
            return(FinalScores)
        }
    }
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    genes_list <- lapply(genes_list, function(g) {
        intersect(genes_list[[g]], rownames(object))
    }) %>% list_clean()
    # 2 extract the data
    data_f <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()

    # 3 run the program
    df_JASMINE <- lapply(seq_along(genes_list), function(g_l) {
        oddsratio <- JASMINE(data_f, genes_list[[g_l]], method = "oddsratio")
        likelihood <- JASMINE(data_f, genes_list[[g_l]], method = "likelihood")
        res_JASMINE <- cbind(oddsratio[, 2], likelihood[, 2])
        colnames(res_JASMINE) <- c(paste0(names(genes_list)[g_l], "_jasmine_or"), paste0(names(genes_list)[g_l], "_jasmine_ll"))
        return(res_JASMINE)
    }) %>%
        do.call(cbind, .) %>%
        as.data.frame()
    object@meta.data <- cbind(object@meta.data, df_JASMINE)
    return(object)
}
