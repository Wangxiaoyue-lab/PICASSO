suppressPackageStartupMessages(require(pacman))
p_load(Seurat, tidyverse)

print(paste("This scripts is in:", getwd()))


ifelse(!exists("work_dir"), warning("Variable 'work_dir' does not exist."),
    print(paste("Your work directory is in:", work_dir))
)
ifelse(!exists("project_name"), warning("Variable 'project_name' does not exist."),
    print(paste("Your project name is:", project_name))
)



script_path <- getwd()
source("./01qc_by_seurat/sc_process.R")



save_file <- function(
    file = NULL, data = NULL, fun = NULL, name_string = NULL, ...) {
    # determine file extension based on write function (if provided)

    defaultend <- switch(deparse(substitute(fun)),
        "saveRDS" = ".rds",
        "write.csv" = ".csv",
        "pdf" = ".pdf"
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
        if (defaultend == ".pdf") {
            fun(filename, ...) # pdf only
        } else {
            return(filename)
        }
    }
}


name_file <- function(work_dir, project_name, name_string = NULL, defaultend = NULL) {
    paste0(work_dir, project_name, "_", name_string, defaultend)
}


find_assay <- function(all_data = all_data) {
    print("The current assay used for the analysis is:", DefaultAssay(all_data))
    print("All assays in this object is:", Assays(all_data))
}


read_refdata <- function(species, file_type) {
    # Check the file_type argument and get the appropriate file path based on the species and file_type arguments
    data <- switch(file_type,
        "cell_cycle_markers" = switch(species,
            "mm" = paste0(script_path, "/Refgenome/cell_cycle_Mus_musculus.csv"),
            "hs" = paste0(script_path, "/Refgenome/cell_cycle_Homo_sapiens.csv"),
            stop("Invalid species argument")
        ),
        "annotations" = switch(species,
            "mm" = paste0(script_path, "/Refgenome/annotations_Mus_musculus.csv"),
            "hs" = paste0(script_path, "/Refgenome/annotations_Homo_sapiens.csv"),
        ),
        stop("Invalid file_type argument")
    ) %>% read.csv()
    return(data)
}

in_data_markers <- function(genes, dataset) {
    not_in_data_marker <- base::setdiff(genes, row.names(dataset))
    if (length(not_in_data_marker) != 0) {
        cat(paste("The following genes are not in the", deparse(substitute(dataset)), ":\n", paste(not_in_data_marker, collapse = ", ")))
        genes <- base::setdiff(genes, not_in_data_marker)
    }
    if (is.null(genes)) {
        stop(paste("No marker genes are in the", deparse(substitute(dataset))))
    }
    return(genes)
}


time_it <- function(f) {
    function(...) {
        start_time <- Sys.time()
        result <- f(...)
        end_time <- Sys.time()
        print(paste0("Execution time: ", end_time - start_time))
        return(result)
    }
}

check_expression <- function(all_data = all_data, feats = NULL, ...) {
    if (any(class(feats) == "data.frame") == TRUE) {
        classed_genes <- names(feats)
        result <- classed_genes %>%
            map(~ {
                mediate_result <- AverageExpression(all_data, features = as.character(na.omit(feats[[.x]]), ...)) %>%
                    as.data.frame() %>%
                    mutate(row_sum = rowSums(.)) %>%
                    rowwise() %>%
                    mutate(row_var = var(c_across(dplyr::everything()))) %>%
                    ungroup()
                mediate_result$class <- .x
                mediate_result
            }) %>%
            bind_rows()
    } else {
        result <- AverageExpression(all_data, features = na.omit(feats), ...) %>%
            as.data.frame() %>%
            mutate(row_sum = rowSums(.)) %>%
            rowwise() %>%
            mutate(row_var = var(c_across(dplyr::everything()))) %>%
            ungroup()
    }
    return(result)
}
