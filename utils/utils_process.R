# 检查列表每个元素长度，如果长于阈值则分裂
list_shorten <- function(object, ...) {
    UseMethod(generic = "list_shorten", object = object)
}
list_shorten.list <- function(object, max_len, name = NULL) {
    n_ <- names(object)
    res <- lapply(seq_along(object), function(obj) {
        list_shorten(object = object[[obj]], max_len = max_len, name = names(object)[obj])
    })
    names(res) <- n_
    return(res)
}
list_shorten.default <- function(object, max_len, name = NULL) {
    object_list <- split(
        object,
        ceiling(seq_along(object) / max_len)
    )
    name_ <- name %||% names(object) %||% deparse(substitute(object))
    names(object_list) <- paste(name_, seq_along(object_list), sep = "_")
    return(object_list)
}


# 去除list和其中元素中的"",NA,NULL,logical(0)或者读取导致的特殊符号
list_clean <- function(object, ...) {
    UseMethod(generic = "list_clean", object = object)
}
list_clean.list <- function(object) {
    lapply(object, function(obj) {
        list_clean(obj)
    }) %>% Filter(Negate(function(x) {
        (is.logical(x[1]) && (length(x) == 0)) || is.null(x) || is.na(x)
    }), .)
}
list_clean.default <- function(object) {
    object <- object[!is.na(object)] 
    object %<>%  gsub(" ","",.) %>%                                
                    gsub("c\\(","",.) %>% 
                        gsub("\\\\","",.) %>% 
                            gsub("\\)","",.) %>% 
                                gsub("\"","",.) %>%
                                    keep(~ nchar(.) >= 1)
    return(object)
}


# list 彻底拉平并保留所有元素名字
list_flat <- function(object) {
    object %<>% flatten
    bool_ <- sapply(object, is.list) %>% any()
    if (bool_ == T) {
        return(list_flat(object))
    } else {
        return(object)
    }
}


# list嵌套拼接名字
list_add_names <- function(object, prefix = "") {
    if (is.list(object)) {
        if (prefix == "") {
            new_prefix <- names(object)
        } else {
            new_prefix <- paste(prefix, names(object), sep = "_")
        }
        names(object) <- new_prefix
        object <- mapply(function(elem, name) {
            list_add_names(elem, prefix = name)
        }, object, new_prefix, SIMPLIFY = FALSE)
    }
    return(object)
}





# 标准化
scale_minmax <- function(...) {
    next
}
