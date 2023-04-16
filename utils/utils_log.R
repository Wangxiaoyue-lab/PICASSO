

if (!exists("work_dir")) {
    warning("Variable 'work_dir' does not exist.")
} else {
    if (!dir.exists(work_dir)) {
        warning("Directory 'work_dir' does not exist.")
        create_dir <- readline(prompt = "Do you want to create a new directory with 'work_dir'? (y/n) ")
        if (tolower(create_dir) == "y") {
            dir.create(work_dir, recursive = T)
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
log_start <- function(...) {
    start_time <<- Sys.time()
    cat("This job starts at:", format(start_time), "\n")
    cat("This script is in:", picasso_path, "\n")
}

# 3 start the log
log_start()

# 4 end the log
log_done <- function(...) {
    end_time <- Sys.time()
    time_difference <- end_time - start_time
    cat("This job ends at:", format(end_time), "\n")
    cat("Total time:", format(time_difference), "\n")
    print(sessionInfo())
}

## when you report a code block
# expr <- expression({...})
log_report <- function(expr, report = T) {
    if (report) {
        start_time <- Sys.time()
        cat("This job starts at:", format(start_time), "\n")
        result <- eval(expr)
        log_done()
        return(result)
    }
}

