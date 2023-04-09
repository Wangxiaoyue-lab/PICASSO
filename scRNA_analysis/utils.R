suppressPackageStartupMessages(require(pacman))
p_load(Seurat, tidyverse)

script_path <- getwd()
cat("This scripts is in:", script_path, "\n")

source("./01_qc_by_seurat/sc_process.R")

if (!exists("work_dir")) {
    warning("Variable 'work_dir' does not exist.")
} else {
    if (!dir.exists(work_dir)) {
        warning("Directory 'work_dir' does not exist.")
        create_dir <- readline(prompt = "Do you want to create a new directory with 'work_dir'? (y/n) ")
        if (tolower(create_dir) == "y") {
            dir.create(work_dir)
            cat("Directory", work_dir, "created.\n")
        }
    }
    cat("Your work directory is in:", work_dir, "\n")
    if (!endsWith(work_dir, "/")) {
        warning("Make sure your directory name ends in a /, otherwise the file name will be wrong.")
    }
}



if (!exists("project_name")) {
    warning("Variable 'project_name' does not exist.")
} else {
    cat("Your project name is:", project_name, "\n")
}



save_file <- function(
    file = NULL, data = NULL, fun = NULL, name_string = NULL, ...) {
    # determine file extension based on write function (if provided)
    defaultend <- switch(deparse(substitute(fun)),
        "saveRDS" = ".rds",
        "write.csv" = ".csv",
        "pdf" = ".pdf"
    )
    name_file <- function(work_dir, project_name, name_string = NULL, defaultend = NULL) {
        paste0(work_dir, project_name, "_", name_string, defaultend)
    }
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







find_assay <- function(all_data = NULL) {
    cat("The current assay used for the analysis is:", DefaultAssay(all_data), "\n")
    cat("All assays in this object:", Assays(all_data), "\n")
}


read_refdata <- function(species, file_type) {
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
        cat("The following genes are not in the", deparse(substitute(dataset)), ":\n", paste(not_in_data_marker, collapse = ", "), "\n")
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
    modi_fun <- function(mediate_result) {
        mediate_result <- as.data.frame(mediate_result)
        original_cols <- colnames(mediate_result)

        mediate_result <- rowwise(mediate_result) %>%
            mutate(
                row_sum = sum(c_across(all_of(original_cols))),
                row_sd = sd(c_across(all_of(original_cols))),
                row_mean = mean(c_across(all_of(original_cols)))
            ) %>%
            ungroup() %>%
            mutate(
                across(all_of(original_cols), list(norm = ~ (.x - row_mean) / row_sd), .names = "norm_{col}")
            )
    }

    if (any(class(feats) == "data.frame") == TRUE) {
        classed_genes <- names(feats)
        result <- classed_genes %>%
            map(~ {
                mediate_result <- AverageExpression(all_data, features = as.character(na.omit(feats[[.x]]), ...)) %>% modi_fun()
                mediate_result$class <- .x
                mediate_result
            }) %>%
            bind_rows()
    } else {
        result <- AverageExpression(all_data, features = na.omit(feats), ...) %>% modi_fun()
    }
    return(result)
}
