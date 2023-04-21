
# start the parallel mode 
utils_parallel <- function(...){
    library(future)
    plan("multisession")
}

# add prefix to cell of seurat object
utils_add_cell_id <- function(object,string=NULL){
    library(rlang)
    string <- string %||% object@meta.data$orig.ident[1]
    new.names <- stringr::str_c(string,colnames(object),sep='_')
    object %<>% RenameCells(new.names=new.names)
    return(object)
}

# merge the list of seurat objects
utils_seurat_merge <- function(seurat_list){
    merge(seurat_list[[1]],seurat_list[-1])
}

utils_read_refdata <- function(species, file_type,picasso_path=picasso_path) {
    data <- switch(file_type,
        "cell_cycle_markers" = switch(species,
            "mm" = paste0(picasso_path,"/other_bioinfo/annotation/Refgenome/cell_cycle_Mus_musculus.csv"),
            "hs" = paste0(picasso_path,"./other_bioinfo/annotation/Refgenome/cell_cycle_Homo_sapiens.csv"),
            stop("Invalid species argument")
        ),
        "annotations" = switch(species,
            "mm" = paste0(picasso_path, "/other_bioinfo/annotation/Refgenome/annotations_Mus_musculus.csv"),
            "hs" = paste0(picasso_path, "/other_bioinfo/annotation/Refgenome/annotations_Homo_sapiens.csv"),
        ),
        stop("Invalid file_type argument")
    ) %>% read.csv()
    return(data)
}