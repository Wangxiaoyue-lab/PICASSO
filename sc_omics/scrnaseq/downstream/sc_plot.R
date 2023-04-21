plot_check_pre <- function(
    object,
    feature_scatter = TRUE,
    vln_group,
    dim_group,
    feats = NULL,
    colors = pal) {
    # @vln_group = "type"
    # @dim_group = c("orig.ident", "type", "doublet_info", "Phase")
    # @feats = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb")
    lapply(vln_group, function(x) {
        print(VlnPlot(object,
            group.by = vln_group,
            features = feats,
            pt.size = 0,
            ncol = 3
        ) + NoLegend())
    })
    message("plot the VlnPlot")
    if (feature_scatter == TRUE) {
        print(QC_Plot_UMIvsGene(object,
            meta_gradient_name = "percent_mito",
            low_cutoff_gene = 1000,
            high_cutoff_UMI = 3000,
            meta_gradient_low_cutoff = 20
        ) +
            scale_x_log10() +
            scale_y_log10())
    }
    message("plot the UMIvsGene scatter plot")
    print(ElbowPlot(object))
    message("plot the ElbowPlot of PCA")
    lapply(dim_group, function(x) {
        print(DimPlot(object, group.by = x))
    })
    message("plot the DimPlot of specified factor")
    lapply(feats, function(x) {
        print(FeaturePlot_scCustom(object, colors_use = colors, features = x))
    })
    message("plot the FeaturePlot of specified continuous variable")
}

plot_processed <- function(
    object,
    dim_group = c("orig.ident", "Phase"), # "type",
    feats = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_hb"),
    colors = pal,
    resolutions = c(0.1, 0.2, 0.3, 0.5)) {
    assay_use <<- ifelse(exists("assay_use"), assay_use %||%
        DefaultAssay(object), DefaultAssay(object))
    lapply(dim_group, function(x) {
        print(DimPlot_scCustom(object, group.by = x, ggplot_default_colors = TRUE, figure_plot = TRUE))
    })
    lapply(feats, function(x) {
        print(FeaturePlot_scCustom(object, colors_use = colors, features = x))
    })

    library(clustree)
    lapply(resolutions, function(x) {
        group <- paste0(assay_use, "_snn_res.", x)
        print(DimPlot(object, group.by = group, label = T))
    })
    print(clustree(object@meta.data, prefix = paste0(assay_use, "_snn_res.")))
}


plot_makers <- function(
    markers,
    object,
    dot_plot = TRUE,
    max_markers = 40,
    cluster = 10,
    feature_plot = TRUE,
    col_dot = viridis_plasma_dark_high,
    col_fea = pal, ...) {
    assay_use <<- ifelse(exists("assay_use"), assay_use %||%
        DefaultAssay(object), DefaultAssay(object))
    plot_function <- function(markers, annotation) {
        # Create a function to plot the dot plot
        draw_dot_plot <- function() {
            # Calculate the number of markers and the number of plots
            n_markers <- length(markers)
            n_plots <- ceiling(n_markers / max_markers)
            # Loop through each plot
            for (i in 1:n_plots) {
                # Calculate the start and end index for the markers
                start <- (i - 1) * max_markers + 1
                end <- min(i * max_markers, n_markers)
                # Plot the dot plot
                if (is.numeric(cluster)) {
                    # Create a vector of markers that are present in the data
                    scaled_markers <- switch(assay_use,
                        "RNA" =  check_markers(markers, object[["RNA"]]@scale.data),
                        "SCT" =  check_markers(markers, object[["SCT"]]@scale.data)
                    )
                    Clustered_DotPlot(object,
                        features = scaled_markers[start:end], k = cluster,
                        plot_km_elbow = FALSE, ...
                    )
                } else {
                    p_dot <- DotPlot_scCustom(object, features = markers[start:end], flip_axes = TRUE, colors_use = col_dot, ...)
                }
                # Add annotation
                if (annotation == TRUE) {
                    print(p_dot + plot_annotation(paste0(cell_type), theme = theme(
                        plot.title = element_text(size = 18, face = "bold")
                    )))
                } else {
                    print(p_dot)
                }
            }
        }

        # Create a function to plot features
        draw_feature_plot <- function() {
            # For each marker
            for (j in markers) {
                # Plot the feature with FeaturePlot_scCustom
                p_fea <- FeaturePlot_scCustom(object,
                    features = j, colors_use = col_fea,
                    ...
                ) & NoAxes()
                # If annotation is true
                if (annotation == TRUE) {
                    # Print the feature plot with the cell type as the title
                    print(p_fea + plot_annotation(paste0(cell_type), theme = theme(
                        plot.title = element_text(size = 18, face = "bold")
                    )))
                } else {
                    # Print the feature plot without annotation
                    print(p_fea)
                }
            }
        }


        # Create the dot plot (if requested)
        if (dot_plot == TRUE) {
            draw_dot_plot()
        }

        # Create the feature plot (if requested)
        if (feature_plot == TRUE) {
            draw_feature_plot()
        }
    }

    # If the unfiltered_markers is a data frame, we have multiple cell types, so we'll loop through them.
    if (class(markers) == "list") {
        # Load the patchwork package so we can annotate plots
        library(patchwork)
        # Get the cell type names from the column names.
        cell_types <- colnames(markers)
        # Loop through the cell types
        for (cell_type in cell_types) {
            # Plot the marker genes for each cell type.
            cat("Ploting markers in ", cell_type, "...\n")
            plot_function(check_markers(markers[, cell_type] %>%
                na.omit() %>% unlist() %>% as.character(), object), annotation = TRUE)
        }
    } else {
        # Otherwise, we have a single cell type, so just plot the marker genes.
        cat("Ploting markers...\n")
        plot_function(check_markers(na.omit(markers), object), annotation = FALSE)
    }
}





plot_celltype_proportions <- function(object, celltype = NULL, group_var = NULL) {
    # Extract the cell type data from the Seurat object
    celltype_data <- object$celltype

    # Calculate the proportions of each cell type
    df <- data.frame(
        clu = names(table(celltype_data)),
        per = sprintf("%1.2f%%", 100 * table(celltype_data) / length(celltype_data))
    )

    # Add the proportion data to the Seurat object
    object$per <- df[match(celltype_data, df$clu), 2]
    object$celltypeper <- paste0(celltype_data, ": (", object$per, ")")

    # Create the stacked bar plot
    if (is.null(group_var)) {
        ggplot(object, aes(x = "", fill = celltypeper)) +
            geom_bar(width = 1) +
            coord_polar("y") +
            theme_void()
    } else {
        ggplot(object, aes(x = group_var, fill = celltypeper)) +
            geom_bar(position = "fill") +
            scale_y_continuous(labels = scales::percent) +
            theme(axis.text.x = element_text(angle = 90))
    }
}
