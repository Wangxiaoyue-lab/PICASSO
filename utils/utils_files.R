# 5 file source system

# save your file including rds,csv and pdf
save_file <- function(file_name = NULL, data = NULL, fun = NULL, sub_dir = NULL, name_string = NULL, ...) {
    # Determine file extension based on write function (if provided)
    default.end <- switch(deparse(substitute(fun)),
        "saveRDS" = ".rds",
        # "save" = ".rdata",
        "write.csv" = ".csv",
        "pdf" = ".pdf"
    )
    # Create subdirectory if it doesn't exist
    if (!is.null(sub_dir)) {
        if (!dir.exists(sub_dir)) {
            warning("'sub_dir' does not exist, creating a new directory.")
            dir.create(paste0(work_dir, "/", sub_dir, "/"), recursive = T)
        }
        if (!endsWith(sub_dir, "/")) {
            sub_dir <- paste0(sub_dir, "/")
        }
    }
    # Use provided file name or generate one
    if (!is.null(file_name)) {
        filename <- file_name
    } else {
        filename <- paste0(work_dir, sub_dir, project_name, "_", name_string, default.end)
    }
    # Write data to file
    if (!is.null(data)) {
        if (default.end == ".pdf") {
            fun(filename, ...) # pdf only
            data %>% lapply(., print)
            dev.off()
        }
        fun(data, filename, ...)
    } else {
        if (default.end == ".pdf") {
            fun(filename, ...) # pdf only
        } else {
            return(filename)
        }
    }
}


read_refdata <- function(species, file_type) {
    data <- switch(file_type,
        "cell_cycle_markers" = switch(species,
            "mm" = paste0(script_path, "../annotation/Refgenome/cell_cycle_Mus_musculus.csv"),
            "hs" = paste0(script_path, "../annotation/Refgenome/cell_cycle_Homo_sapiens.csv"),
            stop("Invalid species argument")
        ),
        "annotations" = switch(species,
            "mm" = paste0(script_path, "../annotation/Refgenome/annotations_Mus_musculus.csv"),
            "hs" = paste0(script_path, "../annotation/Refgenome/annotations_Homo_sapiens.csv"),
        ),
        stop("Invalid file_type argument")
    ) %>% read.csv()
    return(data)
}
