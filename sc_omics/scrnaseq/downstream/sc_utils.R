
utils_add_cell_id <- function(object,string=NULL){
    library(rlang)
    string <- string %||% object@meta.data$orig.ident[1]
    colnames(object) %<>% stringr::str_c(string,.,sep='_')
    return(object)
}
