# report the key common variable
# if (!exists("work_dir")) {
#    warning("Variable 'work_dir' does not exist.")
# } else {
#    if (!dir.exists(work_dir)) {
#        warning("Directory 'work_dir' does not exist.")
#        create_dir <- readline(prompt = "Do you want to create a new directory with 'work_dir'? (y/n) ")
#        if (tolower(create_dir) == "y") {
#            dir.create(work_dir, recursive = T)
#            cat("Directory", work_dir, "created.\n")
#        }
#    }
#    cat("Your work directory is in:", work_dir, "\n")
#    if (!endsWith(work_dir, "/")) {
#        work_dir <<- paste0(work_dir, "/")
#    }
# }

# if (!exists("project_name")) {
#    warning("Variable 'project_name' does not exist.")
# } else {
#    cat("Your project name is:", project_name, "\n")
# }



# report the job
## when you report a long script
log_start <- function(...) {
    start_time <<- Sys.time()
    cat("This job starts at:", format(start_time), "\n")
    cat("The PICASSO path is:", picasso_path, "\n")
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

log_today <- function(...) {
    Sys.Date() %>% format("%Y%m%d")
}


# report the memory
log_memory <- function() {
    library(magrittr)
    mem_list <- system("free -m", intern = T) %>%
        stringr::str_split(pattern = "\\s") %>%
        lapply(function(s) {
            s[nchar(s) > 0]
        })
    mem <- mem_list[[2]][2:4] %>%
        as.numeric() %>%
        (function(x) x / 1024) %>%
        round(., 0)
    names(mem) <- mem_list[[1]][1:3]
    return(mem)
}

# repory the random seed
log_seed <- function(seed = NULL) {
    require(rlang)
    seed <- seed %||% 123 # please don't modify 123!
    log_message(stringr::str_c("The current seed is ", seed))
    set.seed(seed)
}

log_message <- function(..., verbose = T) {
    message(...)
}


log_seurat <- function(object) {
    log_message("The size of seurat object")
    print(format(object.size(object), units = "Mb"))
    log_message("The metadata of seurat object")
    print(str(object@meta.data))
    log_message("The number of seurat object")
    print(length(object@assays))
    print(names(object@assays))
    log_message("The situation of RNA assay")
    print(str(object[["RNA"]]))
    log_message("The situation of reductions")
    print(length(object@reductions))
    print(names(object@reductions))
    log_message("The situation of Commands")
    print(Command(object))
}

log_file <- function(filepath) {
    file_names <- list.files(
        path = filepath,
        # pattern =  NULL,
        recursive = T, full = T
    )
    assertthat::assert_that(!is.null(file_names))
    file_info <- file.info(file_names)[, c("ctime", "size")]
    file_sha256 <- sapply(file_names, function(f) {
        log_sha256(f)
    })
    result <- data.frame(
        name = basename(file_names),
        fullname = file_names,
        create_time = file_info$ctime,
        size = file_info$size,
        sha256 = file_sha256
    )
    print(result)
    return(result)
}


log_sha256 <- function(x) {
    digest::digest(x, file = T, algo = "sha256")
}
