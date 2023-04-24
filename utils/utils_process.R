# 检查列表每个元素长度，如果长于阈值则分裂
list_shorten <- function(object, ...) {
  UseMethod(generic = "list_shorten", object = object)
}
 
list_shorten.list <- function(object,max_len,name=NULL){
    n_ <- names(object)
    res <- lapply(seq_along(object),function(obj){
        list_shorten(object=object[[obj]],max_len=max_len,name=names(object)[obj]) 
    })
    names(res) <- n_
    return(res)
}

list_shorten.default <- function(object,max_len,name=NULL){
    object_list <- split(object,
                    ceiling(seq_along(object)/max_len)) 
    name_ <- name %||% names(object) %||% deparse(substitute(object))                
    names(object_list) <-  paste(name_,seq_along(object_list),sep='_')
    return(object_list)
}


#去除list和其中元素中的"",NA,NULL,logical(0)
list_clean <- function(object, ...) {
  UseMethod(generic = "list_clean", object = object)
}

list_clean.list <- function(object){
    lapply(object,function(obj){
        list_clean(obj)
    }) %>% Filter(Negate(function(x) 
                            (is.logical(x[1]) && (length(x)==0))  || is.null(x) || is.na(x)
                     ),.)
}

list_clean.default <- function(object){
    object <- object[!is.na(object)]
    object %<>% keep(~nchar(.)>=1)
    return(object)
}





scale_minmax <- function(...){
    next
}