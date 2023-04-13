


read_refdata <- function(species,file_type) {
    data <- switch(file_type,
        "cell_cycle_markers" = switch(species,
            "mm" = paste0(script_path, "../annotation/Refgenome/cell_cycle_Mus_musculus.csv"),
            "hs" = paste0(script_path, "../annotation/Refgenome/cell_cycle_Homo_sapiens.csv"),
            stop("Invalid species argument")
        ),
        "annotations" = switch(species,
            "mm" = paste0(script_path, "../annotation/Refgenome/annotations_Mus_musculus.csv"),
            "hs" = paste0(script_path, "../annotation/Refgenome/annotations_Homo_sapiens.csv"),
        ),
        stop("Invalid file_type argument")
    ) %>% read.csv()
    return(data)
}
 