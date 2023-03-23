suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(patchwork))

script_path <- getwd()
# source("sc_process.R")

save_file <- function(
    file = NULL, data = NULL, fun = NULL, name_string = NULL, ...) {
    # determine file extension based on write function (if provided)
    if (!is.null(fun)) {
        switch(deparse(substitute(fun)),
            "saveRDS" = defaultend <- ".rds",
            "write.csv" = defaultend <- ".csv",
            "pdf" = defaultend <- ".pdf",
        )
        # generate file name
        if (!is.null(file)) {
            filename <- file
        } else {
            filename <- name_file(work_dir, project_name, name_string, defaultend)
        }
        # write data
        if (!is.null(data)) {
            fun(data, filename, ...)
        } else {
            return(filename)
        }
    } else {
        # if no write function is provided then just return file name
        warning("No function specified!")
    }
}
name_file <- function(work_dir, project_name, name_string = NULL, defaultend = NULL) {
    paste0(work_dir, project_name, "_", name_string, defaultend)
}

if (exists("run_sctransform")) {
    assay_use <- switch(run_sctransform + 1,
        "RNA",
        "SCT"
    )
} else {
    assay_use <- "RNA"
}
species <- "mm"
cell_cycle_genes_file <- switch(species,
    "mm" = paste0(script_path, "/../Refgenome/cell_cycle_Mus_musculus.csv"),
    "hs" = paste0(script_path, "/../Refgenome/cell_cycle_Homo_sapiens.csv")
)
cell_cycle_genes <- read.csv(cell_cycle_genes_file)

annotations_file <- switch(species,
    "mm" = paste0(script_path, "/../Refgenome/annotations_Mus_musculus.csv"),
    "hs" = paste0(script_path, "/../Refgenome/annotations_Homo_sapiens.csv")
)
annotations <- read.csv(annotations_file)

hb_pattern <- switch(species,
    "mm" = "^Hb[^(p)]",
    "hs" = "^HB[^(P)]"
)


in_data_markers <- function(genes, dataset) {
    not_in_data_marker <- setdiff(genes, row.names(dataset))
    warning(paste("The following genes are not in the", deparse(substitute(dataset), ":", paste(not_in_data_marker, collapse = ", "))))
    marker <- setdiff(genes, not_in_data_marker)
    if (is.null(marker)) {
        stop(paste("No marker genes are in the dataset!", deparse(substitute(dataset))))
    }
    return(marker)
}
