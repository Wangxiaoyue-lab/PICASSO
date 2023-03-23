suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(scCustomize))


# This function checks the quality of the data
qc_check <- function(all_data = all_data) {
    if (is.null(cell_cycle_genes_file) || is.null(annotations_file) || is.null(hb_pattern)) {
        stop("Invalid species argument. Must be 'mm' or 'hs'.")
    }



    # Get gene names for Ensembl IDs for each gene
    cell_cycle_markers <- left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

    # Acquire the S and G2M phase genes
    s_g2m_genes <- cell_cycle_markers %>%
        dplyr::filter(phase %in% c("S", "G2/M")) %>%
        group_by(phase) %>%
        summarise(genes = list(gene_name))

    s_genes <- s_g2m_genes$genes[[which(s_g2m_genes$phase == "S")]]
    g2m_genes <- s_g2m_genes$genes[[which(s_g2m_genes$phase == "G2/M")]]


    all_data_processed <- Add_Mito_Ribo_Seurat(all_data, species = species) %>%
        PercentageFeatureSet(hb_pattern, col.name = "percent_hb") %>%
        CellCycleScoring(g2m.features = g2m_genes, s.features = s_genes) %>%
        NormalizeData() %>%
        ScaleData(features = rownames(all_data)) %>%
        FindVariableFeatures() %>%
        RunPCA(verbose = F, npcs = 20) %>%
        RunUMAP(dims = 1:20) %>%
        Find_doublet()
    return(all_data_processed)


    pdf(save_file(fun = pdf, name_string = "before_qc"))
    print(VlnPlot(all_data_processed,
        group.by = "orig.ident",
        features = feats,
        pt.size = 0,
        ncol = 3
    ) + NoLegend())
    print(QC_Plot_UMIvsGene(all_data_processed,
        meta_gradient_name = "percent_mito",
        low_cutoff_gene = 1000,
        high_cutoff_UMI = 3000,
        meta_gradient_low_cutoff = 20
    ) +
        scale_x_log10() +
        scale_y_log10())
    print(ElbowPlot(all_data_processed))
    print(DimPlot(all_data_processed, split.by = "orig.ident", ncol = 2))
    print(DimPlot(all_data_processed, group.by = gruop_var))
    print(DimPlot(all_data_processed, group.by = "doublet_info"))
    print(DimPlot(all_data_processed, group.by = "Phase"))

    for (i in feats) {
        print(FeaturePlot_scCustom(all_data_processed, colors_use = pal, features = i))
    }
    dev.off()
}

Find_doublet <- function(data) {
    suppressPackageStartupMessages(require(DoubletFinder))
    sweep.res.list <- paramSweep_v3(data, PCs = 1:20, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    p <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>%
        as.character() %>%
        as.numeric()
    nExp_poi <- round(0.05 * ncol(data))
    data <- doubletFinder_v3(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
    c <- grep("pANN_", colnames(data@meta.data))
    data@meta.data <- data@meta.data[, -c]
    return(data)
}


qc_process <- function(all_data,
                       dim_use = 20,
                       resolutions = c(0.1, 0.2, 0.3, 0.5),
                       run_harmony = TRUE,
                       max_mt = 20,
                       max_hb = 5,
                       nFeature_RNA_min = 500,
                       nFeature_RNA_max = 7500,
                       vars_to_regress = c("percent_mito", "S.Score", "G2M.Score"),
                       group_by_vars = "orig.ident") {
    all_data <- subset(all_data,
        subset = percent_mito < max_mt &
            percent_hb < max_hb &
            doublet_info == "Singlet" &
            nFeature_RNA > nFeature_RNA_min &
            nFeature_RNA < nFeature_RNA_max
    )
    ident_value <- paste0(assay_use, "_snn_res.", min(resolutions))

    if (run_sctransform) {
        all_data <- SCTransform(all_data, vars.to.regress = vars_to_regress) %>%
            RunPCA()
    } else {
        all_data <- NormalizeData(all_data) %>%
            FindVariableFeatures() %>%
            ScaleData(vars.to.regress = vars_to_regress) %>%
            RunPCA()
    }

    reduction_use <- switch(run_harmony + 1,
        "pca",
        "harmony"
    )

    if (run_harmony) {
        require(harmony)
        all_data <- all_data %>%
            RunHarmony(group.by.vars = group_by_vars, dims.use = 1:dim_use, assay.use = assay_use)
    }

    all_data <- RunUMAP(all_data, reduction = reduction_use, dims = 1:dim_use) %>%
        FindNeighbors(reduction = reduction_use, dims = 1:dim_use) %>%
        FindClusters(resolution = resolutions) %>%
        SetIdent(value = ident_value)

    return(all_data)
    require(clustree)
    pdf(save_file(fun = pdf, name_string = "after_qc"))

    print(DimPlot(all_data, split.by = "orig.ident", ncol = 2))
    print(DimPlot_scCustom(all_data, group.by = gruop_var, ggplot_default_colors = TRUE, figure_plot = TRUE))
    print(DimPlot_scCustom(all_data, group.by = "Phase", ggplot_default_colors = TRUE, figure_plot = TRUE))
    for (i in feats) {
        print(FeaturePlot_scCustom(all_data, colors_use = pal, features = i))
    }
    for (i in resolutions) {
        group <- paste0(assay_use, "_snn_res.", i)
        print(DimPlot(all_data, group.by = group, label = T))
    }
    print(clustree(all_data@meta.data, prefix = paste0(assay_use, "_snn_res.")))
    dev.off()
}
find_markers <- function(all_data, all = TRUE, ident = NULL, loop_var = NULL, ...) {
    suppressPackageStartupMessages(require(scCustomize))
    modi_fun <- function(marker_genes) {
        marker_genes <- Add_Pct_Diff(marker_genes) %>%
            mutate(type = if_else(avg_log2FC > 0, "up", "down")) %>%
            left_join(
                unique(annotations[, c("gene_name", "description")]),
                by = c("gene" = "gene_name")
            )
        return(marker_genes)
    }
    if (assay_use == "SCT") {
        all_data <- PrepSCTFindMarkers(all_data)
    }
    if (!is.null(ident)) {
        all_data <- SetIdent(all_data, value = ident)
    }
    if (all == TRUE) {
        marker_genes <- FindAllMarkers(all_data, assay = assay_use, ...) %>% modi_fun()
    } else {
        if (!is.null(loop_var)) {
            data_list <- lapply(loop_var, function(x) {
                marker_genes <- FindMarkers(all_data, assay = assay_use, subset.ident = x, ...) %>%
                    rownames_to_column("gene") %>%
                    mutate(cluster = rep(x, nrow(.))) %>%
                    modi_fun()
            })
            marker_genes <- do.call(rbind, data_list)
        } else {
            marker_genes <- FindMarkers(all_data, assay = assay_use, ...) %>%
                rownames_to_column("gene") %>%
                modi_fun()
        }
    }
    return(marker_genes)
}
