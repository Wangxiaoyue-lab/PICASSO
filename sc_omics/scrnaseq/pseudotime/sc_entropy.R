sc_pseudotime_scent <- function(object) {
    library(SCENT)
    data(net13Jun12)
    # 0 check
    assertthat::assert_that(class(object) == "Seurat")
    DefaultAssay(object) <- "RNA"
    species_object <- check_species(object, abbr = T)
    if (species_object == "hs") {
        library(org.Hs.eg.db)
        org <- org.Hs.eg.db
    } else if (species_object == "mm") {
        library(org.Mm.eg.db)
        org <- org.Mm.eg.db
    } else {
        stop("Please check the species of the seurat object")
    }
    # 1 get data
    data_import <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "data"
    ) %>%
        t() %>%
        as.data.frame()
    # 2 convert symbol to entrez id
    ENTREZID <- AnnotationHub::mapIds(org,
        keys = rownames(object),
        column = "ENTREZID",
        keytype = "SYMBOL", #' ENTREZID'
        multiVals = "first"
    )
    row.names(data_import) <- ENTREZID
    # 3 run ccat
    data_import <- log2(data_import + 1)
    ccat.score <- CompCCAT(exp = data_import, ppiA = net13Jun12.m)
    object@meta.data$ccat_score <- ccat.score
    return(object)
}
