# core packages
libraries <- c("scCustomize", "Seurat", "BPCells", "SeuratObject", "SeuratDisk")
lapply(
    libraries,
    function(x) {
        suppressPackageStartupMessages(library(x, character.only = TRUE))
    }
)



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
            "mm" = paste0(picasso_path, "/other_bioinfo/annotation/reference/cell_cycle_Mus_musculus.csv"),
            "hs" = paste0(picasso_path, "./other_bioinfo/annotation/reference/cell_cycle_Homo_sapiens.csv"),
            stop("Invalid species argument")
        ),
        "annotations" = switch(species,
            "mm" = paste0(picasso_path, "/other_bioinfo/annotation/reference/annotations_Mus_musculus.csv"),
            "hs" = paste0(picasso_path, "/other_bioinfo/annotation/reference/annotations_Homo_sapiens.csv"),
        ),
        stop("Invalid file_type argument")
    ) %>% read.csv()
    return(data)
}

utils_markers_scb <- function(path = NULL, species, clean = F, disease_pattern = NULL) {
    path <- path %||% "/other_bioinfo/annotation/annotation_files/marker"
    assertthat::assert_that(species %in% c("hs", "mm"))
    path_s <- switch(species,
        "hs" = paste0(picasso_path, path, "/singleCellBase_human_cell_markers_20220629.txt"),
        "mm" = paste0(picasso_path, path, "/singleCellBase_mouse_cell_markers_20220629.txt")
    )
    marker_df <- data.table::fread(path_s) %>%
        as.data.frame()
    if (clean == T) {
        marker_df %<>%
            as_tibble() %>%
            select(
                cell_type, gene_symbol, sample_type,
                disease, sci_if, date
            ) %>%
            filter(grepl(disease, pattern = disease_pattern)) %>%
            select(cell_type, gene_symbol) %>%
            arrange(cell_type) %>%
            group_by(cell_type) %>%
            nest() %>%
            mutate(markers = map(data, ~ {
                tryCatch(
                    {
                        paste0(.x, collapse = ",") %>%
                            str_split(., pattern = ",", simplify = T) %>%
                            .[1, ] %>%
                            unique()
                    },
                    error = function(e) {
                        NA
                    }
                )
            })) %>%
            select(cell_type, markers) %>%
            deframe() %>%
            list_clean()
    }
    return(marker_df)
}
