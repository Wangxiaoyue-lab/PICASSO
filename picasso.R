libraries <- c("dplyr")
lapply(
    libraries,
    function(x) {
        suppressPackageStartupMessages(library(x, character.only = TRUE))
    }
)

picasso_path <- getwd()

# file management system functions: choose_pipeline, list_pipeline, load_necessary and load_script

# choose the pipeline
choose_pipeline <- function(pipeline = NULL,
                            module = NULL) {
    # If the pipeline name is not specified, list available pipelines
    if (is.null(pipeline)) {
        list_pipeline()
        pipeline <- readline(prompt = "Please specify a pipeline and press 'Enter'\n")
    }

    # Load all necessary packages and functions
    load_necessary()

    # Load scripts for the selected pipeline
    lapply(pipeline, function(p) {
        list.files(path = picasso_path, recursive = F, full = T) %>%
            lapply(., function(ls) {
                pipe <- list.files(path = ls, recursive = F, full = T) %>%
                    grep(., pattern = p, value = T)
                if (!is.null(module)) {
                    pipe <- lapply(module, function(m) {
                        list.dirs(path = pipe, recursive = F, full = T) %>%
                            grep(., pattern = m, value = T)
                    })
                }
            })
    }) %>%
        unlist(., recursive = T) %>%
        stringr::str_split(., pattern = "PICASSO/", simplify = T, n = 2) %>%
        .[, 2] %>%
        lapply(., function(rs) {
            rs %>% load_script(dir = .)
        })

    # Print a message indicating the pipeline has been loaded
    cat(paste("Succeed to load script:", pipeline, collapse = "\n"), "\n")
}

# list all pipelines
list_pipeline <- function(pipeline = NULL, module = F) {
    # Get all of the pipelines
    dir_1 <- list.dirs(path = picasso_path, recursive = F)
    black_list <- c("\\.git", "\\.vscode", "picture", "knowledge_base", "utils")
    dir_1 <- dir_1[!grepl(dir_1, pattern = paste(black_list, collapse = "|"))]
    for (d in dir_1) {
        # Get the name of the total class
        total_class <- stringr::str_split(d, pattern = "PICASSO/", simplify = T, n = 2)[, 2]
        cat("\n#----", total_class, "----#\n")
        # Get all of the pipelines within the total class
        dir_2 <- list.dirs(path = d, recursive = F)
        for (p in dir_2) { # pipeline
            # Get the name of the pipeline
            pipe_exist <- stringr::str_split(p, pattern = paste0(total_class, "/"), simplify = T, n = 2)[, 2]
            if (is.null(pipe_exist)) {
                next
            }
            # If the user specifies a pipeline, only include it in the list
            if (!is.null(pipeline)) {
                pipe_exist_bool <- grepl(pipe_exist, pattern = pipeline)
                if (!pipe_exist_bool) {
                    next
                }
            }
            cat("-->", pipe_exist, "\n")
            # Get all of the modules within the pipeline
            dir_3 <- list.dirs(path = grep(p, pattern = pipe_exist, value = T), recursive = F)
            if (module == T) {
                for (m in dir_3) { # module
                    # Get the name of the module
                    modules <- stringr::str_split(m, pattern = paste0(pipe_exist, "/"), simplify = T, n = 2)[, 2]
                    if (is.null(modules)) {
                        next
                    }
                    cat("--> -->", modules, "\n")
                }
            }
        }
    }
}



# load necessary
load_necessary <- function(...) {
    ## basical utils
    load_script(dir = "utils", script = "core_utils")
    load_script(dir = "utils/parallel", script = "parallel")
    # load_script(dir = "utils/input_your_parameter", script = "parameter")
    ## color
    # load_script(dir='visualization/colour',script='palette')
    ## picture
    # load_script(dir='visualization/plot',script='themes')
}


# load the specified script
load_script <- function(dir, script = "\\.R") {
    list.files(
        path = paste0(picasso_path, "/", dir),
        pattern = paste(script, collapse = "|"),
        recursive = T, full = T
    ) %>%
        grep(., pattern = "\\.R", value = T) %>%
        lapply(., source)
}
