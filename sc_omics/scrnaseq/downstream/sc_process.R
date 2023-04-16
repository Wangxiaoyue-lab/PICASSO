process_process <- function(object,
                            npcs = 20,
                            resolutions = c(0.1, 0.2, 0.3, 0.5),
                            future = FALSE,
                            run_harmony = TRUE,
                            run_sctransform = TRUE,
                            group_in_harmony = "orig.ident",
                            vars_to_regress = c("percent_mito", "S.Score", "G2M.Score")) {
    if (future) {
        options(future.globals.maxSize = 1e9)
        plan("multisession")
    }
    assay_use <<- switch(run_sctransform + 1,
        "RNA",
        "SCT"
    )
    ident_value <- paste0(assay_use, "_snn_res.", min(resolutions))
    if (run_sctransform == TRUE) {
        object <- SCTransform(object, vars.to.regress = vars_to_regress) %>%
            RunPCA()
    } else {
        object <- NormalizeData(object) %>%
            FindVariableFeatures() %>%
            ScaleData(features = row.names(object), vars.to.regress = vars_to_regress) %>%
            RunPCA()
    }

    reduction_use <- switch(run_harmony + 1,
        "pca",
        "harmony"
    )

    if (run_harmony) {
        p_load(harmony)
        object <- object %>%
            RunHarmony(group.by.vars = group_in_harmony, dims.use = 1:npcs, assay.use = assay_use)
    }
    object <- RunUMAP(object, reduction = reduction_use, dims = 1:npcs) %>%
        FindNeighbors(reduction = reduction_use, dims = 1:npcs) %>%
        FindClusters(resolution = resolutions) %>%
        SetIdent(value = ident_value)

    return(object)
}


# process_integration



process_find_markers <- function(object,
                                 is_all = TRUE,
                                 future = FALSE,
                                 ident = NULL,
                                 loop_var = NULL,
                                 species, ...) {
    # @is_all：FindAllMarkers or FindMarkers
    # @loop_var: compare conditions within some idents

    if (future) {
        options(future.globals.maxSize = 1e9)
        plan("multicore")
    }
    assay_use <<- ifelse(exists("assay_use"), assay_use %||%
        DefaultAssay(object), DefaultAssay(object))

    add_fun <- function(marker_genes) {
        annotations <- read_refdata(species, "annotations")
        # Add a column with the percentage difference
        marker_genes <- Add_Pct_Diff(marker_genes) %>%
            # Create a new column with the type of marker gene
            mutate(type = if_else(avg_log2FC > 0, "up", "down")) %>%
            # Use the left_join() function to add the gene name and description
            left_join(
                unique(annotations[, c("gene_name", "description")]),
                by = c("gene" = "gene_name")
            )
        return(marker_genes)
    }
    if (assay_use == "SCT") {
        object <- PrepSCTFindMarkers(object)
    }
    if (!is.null(ident)) {
        object <- SetIdent(object, value = ident)
    }
    if (is_all == TRUE) {
        # If all is TRUE, find all markers
        marker_genes <- FindAllMarkers(object, assay = assay_use, ...) %>% add_fun()
    } else {
        if (!is.null(loop_var)) {
            # If all is FALSE and loop_var is not NULL, find markers for each cluster
            data_list <- lapply(loop_var, function(x) {
                marker_genes <- FindMarkers(object, assay = assay_use, subset.ident = x, ...) %>%
                    rownames_to_column("gene") %>%
                    mutate(cluster = rep(x, nrow(.))) %>%
                    add_fun()
            })
            marker_genes <- do.call(rbind, data_list)
        } else {
            # If all is FALSE and loop_var is NULL, find markers
            marker_genes <- FindMarkers(object, assay = assay_use, ...) %>%
                rownames_to_column("gene") %>%
                add_fun()
        }
    }
    return(marker_genes)
}


# process_celltype


# process_individual_deg