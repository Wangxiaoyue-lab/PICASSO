suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(patchwork))

source("sc_process.R")

save_files <- function(
    file = NULL, data = NULL, fun = NULL, name_string = NULL, ...) {
    # determine file extension based on write function (if provided)
    if (!is.null(fun)) {
        switch(deparse(substitute(fun)),
            "saveRDS" = defaultend <- ".rds",
            "write.csv" = defaultend <- ".csv",
            "pdf" = defaultend <- ".pdf",
        )

        # generate final file name
        if (!is.null(file)) {
            filename <- file
        } else {
            filename <- name_file(work_dir, project_name, name_string, defaultend)
        }

        # write data
        if (!is.null(data)) {
            fun(data, filename, ...)
        }
    } else {
        # if no write function is provided then just return final file name
        name_file(work_dir, project_name, name_string, defaultend)
        warning("No function specified!")
    }
    name_file <- function(work_dir, project_name, name_strings = name_string, defaultend = NULL) {
        paste0(work_dir, project_name, "_", name_string, defaultend)
    }
}
