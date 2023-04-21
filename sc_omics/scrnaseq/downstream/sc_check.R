# core packages
libraries <- c("scCustomize", "Seurat")
lapply(
    libraries,
    function(x) {
        suppressPackageStartupMessages(library(x, character.only = TRUE))
    }
)

# source("../utils/load_ref.R")


## check the assays
check_assay <- function(object = NULL) {
    cat("The current assay used for the analysis is:", DefaultAssay(object), "\n")
    cat("All assays in this object:", Assays(object), "\n")
}



## check whether the markes exists in the features of object
check_markers <- function(genes, object) {
    not_in_data_marker <- base::setdiff(genes, row.names(object))
    if (length(not_in_data_marker) != 0) {
        cat("The following genes are not in the", deparse(substitute(object)), ":\n", paste(not_in_data_marker, collapse = ", "), "\n")
        genes <- base::setdiff(genes, not_in_data_marker)
    }
    if (is.null(genes)) {
        stop(paste("No marker genes are in the", deparse(substitute(object))))
    }
    return(genes)
}


# This function checks the quality of the data
check_pre <- function(
    object,
    species = c("hs", "mm"),
    cell_cycle_source = c("seurat", "local"),
    npcs = 20,
    # check_doublet = TRUE,
    verbose = F) {
    hb_pattern <- switch(species,
        "mm" = "^Hb[^(p)]",
        "hs" = "^HB[^(P)]",
    )
    if (cell_cycle_source == "seurat") {
        s_genes <- switch(species,
            "hs" = cc.genes$s.genes,
            "mm" = cc.genes$s.genes %>% str_to_tittle()
        )
        g2m_genes <- switch(species,
            "hs" = cc.genes$g2m.genes,
            "mm" = cc.genes$g2m.genes %>% str_to_tittle()
        )
    } else {
        # Get gene names for Ensembl IDs for each gene
        cell_cycle_markers <- left_join(utils_read_refdata(species, "cell_cycle_markers"),
            utils_read_refdata(species, "annotations"),
            by = c("geneID" = "gene_id")
        )
        # Acquire the S and G2M phase genes
        s_g2m_genes <- cell_cycle_markers %>%
            dplyr::filter(phase %in% c("S", "G2/M")) %>%
            group_by(phase) %>%
            summarise(genes = list(gene_name))
        s_genes <- s_g2m_genes$genes[[which(s_g2m_genes$phase == "S")]]
        g2m_genes <- s_g2m_genes$genes[[which(s_g2m_genes$phase == "G2/M")]]
    }
    object <- Add_Mito_Ribo_Seurat(object, species = species) %>%
        PercentageFeatureSet(hb_pattern, col.name = "percent_hb") %>%
        CellCycleScoring(g2m.features = g2m_genes, s.features = s_genes) %>%
        NormalizeData(verbose = verbose) %>%
        ScaleData(features = rownames(object), verbose = verbose) %>%
        FindVariableFeatures(verbose = verbose) %>%
        RunPCA(verbose = verbose, npcs = npcs) %>%
        RunUMAP(dims = 1:npcs, verbose = verbose)
    # if (check_doublet) {
    #    object <- check_doublet(object, npcs)
    # }
    return(object)
}


# check whether doublets exist
check_doublet <- function(object,
                          npcs, # the number of npcs of pca
                          celltype = NULL, # the colname of celltype in metadata
                          ncelltype = NULL, # the accessed number of cell types
                          fast = FALSE # whether use the fast mode(scds package)
) {
    if (fast == F) {
        library(DoubletFinder)
        process_ <- Command(object) %>%
            grepl(., pattern = "PCA") %>%
            Reduce("+", .) %||% 0
        if (process_ < 1) {
            object %<>% NormalizeData(verbose = F) %>%
                ScaleData(features = rownames(object), verbose = F) %>%
                FindVariableFeatures(verbose = F) %>%
                RunPCA(verbose = F, npcs = npcs)
        }
        p <- object %>%
            paramSweep_v3(., PCs = 1:npcs, sct = FALSE) %>%
            summarizeSweep(., GT = FALSE) %>%
            find.pK() %>%
            filter(BCmetric == max(BCmetric)) %>%
            pull(pK) %>%
            as.character() %>%
            as.numeric()
        nExp_poi <- round(ncol(object) * ncol(object) * 1.6 * 1.6 / 2e5)
        if (!is.null(celltype)) {
            homotypic.prop <- object@meta.data %>%
                select(!!sym(celltype)) %>%
                modelHomotypic()
        } else {
            ncelltype <- ncelltype %||% 5
            homotypic.prop <- 1 / ncelltype
        }
        nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
        object <- doubletFinder_v3(object,
            PCs = 1:npcs, pN = 0.25, pK = p, nExp = nExp_poi,
            reuse.pANN = FALSE, sct = FALSE
        )
        colnames(object@meta.data)[ncol(object@meta.data)] <- "doublet_info"
        c <- grep("pANN_", colnames(object@meta.data))
        object@meta.data <- object@meta.data[, -c]
    } else {
        library(scds)
        library(SingleCellExperiment)
        scds_res <- as.SingleCellExperiment(object) %>%
            cxds(., estNdbl = TRUE) %>%
            bcds(., estNdbl = TRUE) %>%
            cxds_bcds_hybrid(., estNdbl = TRUE) %>%
            colData() %>%
            .@listData %>%
            as.data.frame() %>%
            select(cxds_score, cxds_call, bcds_score, bcds_call, hybrid_score, hybrid_call)
        object <- AddMetaData(object, scds_res)
    }
    return(object)
}
# check_doublet <- function(object, npcs) {
#    library(DoubletFinder)
#    sweep.res.list <- paramSweep_v3(object, PCs = 1:npcs, sct = FALSE)
#    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
#    bcmvn <- find.pK(sweep.stats)
#    p <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>%
#        as.character() %>%
#        as.numeric()
#    nExp_poi <- round(0.05 * ncol(object))
#    object <- doubletFinder_v3(object,
#        PCs = 1:npcs, pN = 0.25, pK = p, nExp = nExp_poi,
#        reuse.pANN = FALSE, sct = FALSE
#    )
#    colnames(object@meta.data)[ncol(object@meta.data)] <- "doublet_info"
#    c <- grep("pANN_", colnames(object@meta.data))
#    object@meta.data <- object@meta.data[, -c]
#    return(object)
# }

check_size_future <- function(object) {
    maxSize <- ifelse(dim(object)[2] < 5e4, 1e9,
        ifelse(dim(object)[2] < 1e5, 1e10, 9e10)
    )
    return(maxSize)
}


check_seurat <- function(object) {
    cat("The size is ", dim(object), "\n\n")
    cat("The meta.data is ")
    print(object@meta.data[1:4, ])
    cat("\n\nThe structure is ")
    print(str(object))
}
