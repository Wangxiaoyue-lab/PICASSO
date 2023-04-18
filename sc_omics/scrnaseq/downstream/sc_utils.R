
utils_add_cell_id <- function(object,string=NULL){
    library(rlang)
    string <- string %||% object@meta.data$orig.ident[1]
    new.names <- stringr::str_c(string,colnames(object),sep='_')
    object %<>% RenameCells(new.names=new.names)
    return(object)
}
