#' Check the current and available assays in a Seurat object
#'
#' This function takes a Seurat object and prints the current assay used for analysis as well as all available assays in the object.
#'
#' @param object A Seurat object.
#' @return NULL
#' @examples
#' library(Seurat)
#' pbmc_small <- pbmc_small
#' check_assay(pbmc_small)
check_assay <- function(object = NULL) {
    cat("The current assay used for the analysis is:", DefaultAssay(object), "\n")
    cat("All assays in this object:", Assays(object), "\n")
}



#' Check if marker genes are present in a data object
#'
#' This function takes a character vector of marker genes and a data object (such as a Seurat object or data frame) and checks if the marker genes are present in the row names of the data object. The function prints any marker genes that are not present in the data object and returns a character vector of marker genes that are present in the data object.
#'
#' @param genes A character vector of marker genes.
#' @param object A data object with row names (such as a Seurat object or data frame).
#' @return A character vector of marker genes that are present in the data object.
#' @examples
#' library(Seurat)
#' pbmc_small <- pbmc_small
#' check_markers(c("MS4A1", "CD3E", "NOTAGENE"), pbmc_small)
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


#' Preprocess a Seurat object
#'
#' This function takes a Seurat object and performs several preprocessing steps, including adding mitochondrial and ribosomal gene expression information, cell cycle scoring, normalization, scaling, variable feature selection, PCA and UMAP. The function allows for customization of several parameters.
#'
#' @param object A Seurat object.
#' @param species A character string specifying the species for the analysis ("hs" for human or "mm" for mouse).
#' @param cell_cycle_source A character string specifying the source of cell cycle genes ("seurat" or "local").
#' @param npcs An integer specifying the number of principal components to use for downstream analysis.
#' @param verbose A logical value indicating whether to print progress messages.
#' @return A preprocessed Seurat object.
#' @examples
#' library(Seurat)
#' pbmc_small <- pbmc_small
#' check_pre(pbmc_small)
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



#' Check for doublets in a Seurat object
#'
#' This function takes a Seurat object and performs doublet detection using either the DoubletFinder or scds package. The function allows for customization of several parameters.
#'
#' @param object A Seurat object.
#' @param npcs An integer specifying the number of principal components to use for downstream analysis.
#' @param celltype An optional character string specifying the column name of cell type information in the metadata.
#' @param ncelltype An optional integer specifying the number of cell types to use for homotypic proportion estimation.
#' @param fast A logical value indicating whether to use the fast mode (scds package) for doublet detection.
#' @return A Seurat object with added doublet information in the metadata.
#' @examples
#' library(Seurat)
#' pbmc_small <- pbmc_small
#' check_doublet(pbmc_small, npcs = 10)
check_doublet <- function(object,
                          npcs,
                          celltype = NULL,
                          ncelltype = NULL,
                          fast = FALSE) {
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


#' Check the maximum size of an object
#'
#' @param object The object to check
#' @return The maximum size of the object
#' @export
check_size_future <- function(object) {
    maxSize <- ncol(object) * 4e5
    return(maxSize)
}
