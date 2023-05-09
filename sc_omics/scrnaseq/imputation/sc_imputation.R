# https://github.com/Winnie09/imputationBenchmark
sc_imputation <- function(...) {
    sc_imputation_magic <- function() {
        next
    }
    sc_imputation_knn_smooth <- function() {
        next
    }
    sc_imputation_saver <- function() {
        next
    }
}


sc_imputation_saver <- function(object, ...) {
    UseMethod(generic = "sc_imputation_saver", object = object)
}

sc_imputation_saver.Seurat <- function(object, genes = NULL, n_cores, ...) {
    expression <- as.matrix(object[["RNA"]]@counts)
    result <- sc_imputation_saver(object = expression, genes = genes, n_cores = n_cores, ...)
    return(result)
}

sc_imputation_saver.default <- function(object, genes, n_cores, ...) {
    if (!is.null(genes)) {
        genes <- which(rownames(object) %in% genes)
    }
    saver_genes <- saver(object,
        pred.genes = genes,
        ncores = n_cores,
        estimates.only = TRUE, ...
    )
    return(saver_genes)
}
