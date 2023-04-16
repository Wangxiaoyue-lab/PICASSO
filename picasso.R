require(dplyr)
require(rlang, quietly = T)
picasso_path <- getwd()

choose_pipeline <- function(pipeline = NULL,
                            module = NULL) {
    if (is.null(pipeline)) {
        list_pipeline()
        stop("Please specify a pipeline.")
    }
    load_necessary()
    lapply(pipeline, function(p) {
        list.files(path = picasso_path, recursive = F) %>%
            lapply(., function(ls) {
                pipe <- list.files(path = ls, recursive = F) %>% grep(., pattern = p, value = T)
                if (is.null(module)) {
                    pipe <- pipe
                } else {
                    pipe <- lapply(module, function(m) {
                        list.files(path = paste0(pipe, "/", m), recursive = F)
                    }) %>% unlist()
                }
                return(pipe)
            }) %>%
            unlist()
    }) %>%
        unlist() %>%
        stringr::str_split(., pattern = "PICASSO/", simplify = T, n = 2) %>%
        lapply(., function(rs) {
            load_script(dir = rs, script = "\\.R")
        })
}



list_pipeline <- function(pipeline = NULL, module = F) {
    dir_1 <- list.dirs(path = picasso_path, recursive = F)
    black_list <- c("\\.git", "\\.vscode", "picture", "knowledge_base", "utils")
    dir_1 <- dir_1[!grepl(dir_1, pattern = paste(black_list, collapse = "|"))]
    for (d in dir_1) { # total class
        total_class <- stringr::str_split(d, pattern = "PICASSO/", simplify = T, n = 2)[, 2]
        cat("#----", total_class, "----#\n")
        dir_2 <- list.dirs(path = d, recursive = F)
        for (p in dir_2) { # pipeline
            pipe_exist <- stringr::str_split(p, pattern = paste0(total_class, "/"), simplify = T, n = 2)[, 2]
            if (is.null(pipeline)) {
                pipe_exist <- pipe_exist
            } else {
                pipe_exist <- grep(pipe_exist, pattern = pipeline, value = T)
            }
            cat("-->", pipe_exist, "\n")
            dir_3 <- list.dirs(path = grep(p, pattern = pipe_exist, value = T), recursive = F)
            if (module == T) {
                for (m in dir_3) { # module

                    modules <- stringr::str_split(m, pattern = paste0(pipe_exist, "/"), simplify = T, n = 2)[, 2]
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
    load_script(dir = "utils/input_your_parameter", script = "parameter")
    ## color
    # load_script(dir='visualization/colour',script='palette')
    ## picture
    # load_script(dir='visualization/plot',script='themes')
}


# load the specified script
load_script <- function(dir, script) {
    scripts <- list.files(
        path = paste0(picasso_path, "/", dir),
        pattern = paste(script, collapse = "|"),
        recursive = T, full = T
    )
    lapply(scripts, source)
}
