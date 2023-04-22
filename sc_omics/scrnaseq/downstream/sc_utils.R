#' Start parallel mode
#'
#' @param ... Additional arguments
#'
#' @return NULL
#'
#' @export
utils_parallel <- function(...) {
    library(future)
    plan("multisession")
}

#' Add prefix to cell IDs of a Seurat object
#'
#' @param object A Seurat object
#' @param string A character vector specifying the prefix to add (default: NULL)
#'
#' @return A Seurat object with modified cell IDs
#'
#' @export
utils_add_cell_id <- function(object, string = NULL) {
    library(rlang)
    string <- string %||% object@meta.data$orig.ident[1]
    new.names <- stringr::str_c(string, colnames(object), sep = "_")
    object %<>% RenameCells(new.names = new.names)
    return(object)
}


#' Merge a list of Seurat objects
#'
#' @param seurat_list A list of Seurat objects
#'
#' @return A merged Seurat object
#'
#' @export
utils_seurat_merge <- function(seurat_list) {
    merge(seurat_list[[1]], seurat_list[-1])
}

#' Read reference data
#'
#' @param species A character vector specifying the species
#' @param file_type A character vector specifying the type of file to read
#' @param picasso_path A character vector specifying the path to the Picasso directory (default: picasso_path)
#'
#' @return A data frame containing the reference data
#'
#' @export
utils_read_refdata <- function(species, file_type, picasso_path = picasso_path) {
    data <- switch(file_type,
        "cell_cycle_markers" = switch(species,
            "mm" = paste0(picasso_path, "/other_bioinfo/annotation/Refgenome/cell_cycle_Mus_musculus.csv"),
            "hs" = paste0(picasso_path, "./other_bioinfo/annotation/Refgenome/cell_cycle_Homo_sapiens.csv"),
            stop("Invalid species argument")
        ),
        "annotations" = switch(species,
            "mm" = paste0(picasso_path, "/other_bioinfo/annotation/Refgenome/annotations_Mus_musculus.csv"),
            "hs" = paste0(picasso_path, "/other_bioinfo/annotation/Refgenome/annotations_Homo_sapiens.csv"),
        ),
        stop("Invalid file_type argument")
    ) %>% read.csv()
    return(data)
}
