# 1 core packages
libraries <- c("tidyverse", "data.table","rlang")
lapply(libraries, 
    function(x){
        suppressPackageStartupMessages(library(x, character.only = TRUE))
    }
)

# 2 set the basical variable
if (!exists("work_dir")) {
    warning("Variable 'work_dir' does not exist.")
} else {
    if (!dir.exists(work_dir)) {
        warning("Directory 'work_dir' does not exist.")
        create_dir <- readline(prompt = "Do you want to create a new directory with 'work_dir'? (y/n) ")
        if (tolower(create_dir) == "y") {
            dir.create(work_dir,recursive = T)
            cat("Directory", work_dir, "created.\n")
        }
    }
    cat("Your work directory is in:", work_dir, "\n")
    if (!endsWith(work_dir, "/")) {
        work_dir <<- paste0(work_dir, "/")
    }
}

if (!exists("project_name")) {
    warning("Variable 'project_name' does not exist.")
} else {
    cat("Your project name is:", project_name, "\n")
}



# report the job
## when you report a long script
log_start <- function(...){
    start_time <<- Sys.time()
    cat("This job starts at:", format(start_time), "\n")
    script_path <<- getwd()
    cat("This script is in:", script_path, "\n")
}

# 3 start the log
log_start()

# 4 end the log 
log_done <- function(...){
    end_time <- Sys.time()
    time_difference <- end_time - start_time
    cat("This job ends at:", format(end_time), "\n")
    cat("Total time:", format(time_difference), "\n")
    print(sessionInfo())
}

## when you report a code block
#expr <- expression({...})
log_report <- function(expr,report=T) {
    if(report){
        start_time <- Sys.time()
        cat("This job starts at:", format(start_time), "\n") 
        result <- eval(expr)
        log_done()
        return(result)
    }

}
 

# 5 file source system

# save your file including rds,csv and pdf
save_file <- function(file_name = NULL, data = NULL, fun = NULL, sub_dir = NULL, name_string = NULL, ...) {
    # Determine file extension based on write function (if provided)
    default.end <- switch(deparse(substitute(fun)),
        "saveRDS" = ".rds",
        #"save" = ".rdata",
        "write.csv" = ".csv",
        "pdf" = ".pdf"
    )
    # Create subdirectory if it doesn't exist
    if (!is.null(sub_dir)) {
        if (!dir.exists(sub_dir)) {
            warning("'sub_dir' does not exist, creating a new directory.")
            dir.create(paste0(work_dir, "/", sub_dir, "/"),recursive=T)
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
            data %>% lapply(.,print)
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


