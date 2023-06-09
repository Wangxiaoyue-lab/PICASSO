#' Plot various graphs for checking data quality
#'
#' @param object The object to plot
#' @param feature_scatter Logical, whether to plot UMI vs Gene scatter plot
#' @param vln_group Grouping variable for violin plot
#' @param dim_group Grouping variable for dimensionality reduction plot. i.e. c("orig.ident", "type", "doublet_info", "Phase")
#' @param feats Features to plot in violin and feature plots. i.e. c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb")
#' @param colors Color palette to use
#'
#' @return Various plots for checking data quality
plot_check_pre <- function(
    object,
    feature_scatter = TRUE,
    vln_group,
    dim_group,
    feats = NULL,
    colors = pal) {
    lapply(vln_group, function(x) {
        ncol_ <- object@meta.data %>%
            pull(any_of(vln_group)) %>%
            unique() %>%
            length()
        ncol <- ifelse(ncol_ > 20, 1,
            ifelse(ncol_ > 10, 2, 3)
        )
        print(VlnPlot(object,
            group.by = vln_group,
            features = feats,
            pt.size = 0,
            ncol = ncol
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

plot_features_group <- function(object, features, group, select, pal) {
    library(plot1cell)
    complex_featureplot(object,
        features = features,
        group = group,
        select = select,
        order = T
    )
}

# joint features scatterplot
plot_features_density <- function(object, features, repair = T) {
    # if(repair==T){ #repair the missing gene expression
    Nebulosa::plot_density(object, features, joint = T)
    # }else{
    #    scCustomize::Plot_Density_Joint_Only(seurat_object = object, features = features)
    # }
}

plot_stack_violin <- function(object, features, pal) {
    scCustomize::Stacked_VlnPlot(
        seurat_object = object,
        features = features,
        x_lab_rotate = TRUE,
        colors_use = pal
    )
}

plot_big_circlize <- function(object, core_col, other_col) {
    library(plot1cell)
    prepare <- prepare_circlize_data(object, scale = 0.8)
    cols <- c(core_col, other_col)
    cols_select <- lapply(cols, function(x) {
        object@meta.data %>%
            pull(any_of(x)) %>%
            unique() %>%
            length() %>%
            rand_color()
    }) %>% set_names(cols)
    plot_circlize(prepare,
        do.label = T,
        pt.size = 0.01,
        col.use = cols_select[[core_col]],
        bg.color = "white",
        kde2d.n = 200,
        repel = T,
        label.cex = 0.6
    )
    for (i in 2:length(cols)) {
        add_track(prepare, group = cols[i], colors = cols_select[i], track_num = i)
    }
}


#' Plot various graphs for processed data
#'
#' @param object The object to plot
#' @param dim_group Grouping variable for dimensionality reduction plot
#' @param feats Features to plot in feature plots
#' @param colors Color palette to use
#' @param resolutions Resolutions for clustering
#'
#' @return Various plots for processed data
plot_processed <- function(
    object,
    dim_group = c("orig.ident", "Phase"), # "type",
    feats = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_hb"),
    colors = pal,
    resolutions = c(0.1, 0.2, 0.3, 0.5)) {
    assertthat::assert_that(class(object)=="Seurat")
    assertthat::assert_that(all(dim_group %in% colnames(object@meta.data)))
    assertthat::assert_that(all(feats %in% colnames(object@meta.data)))
    assay_use <<- ifelse(exists("assay_use"), assay_use %||%
        DefaultAssay(object), DefaultAssay(object))
    lapply(dim_group, function(x) {
        print(DimPlot_scCustom(object, group.by = x, ggplot_default_colors = TRUE, figure_plot = TRUE))
    })
    lapply(feats, function(x) {
        print(FeaturePlot_scCustom(object, colors_use = colors, features = x))
    })
    lapply(resolutions, function(x) {
        group <- paste0(assay_use, "_snn_res.", x)
        print(DimPlot(object, group.by = group, label = T))
    })
    library(clustree)
    print(clustree(object@meta.data, prefix = paste0(assay_use, "_snn_res.")))
}


# Plot marker genes
# Suggested plot size: width = 35, height = 20
plot_markers <- function(
    object,
    markers,
    ident_group = NULL,
    dot_plot = TRUE,
    dot_max = 40,
    dot_cluster = NULL,
    dot_col = viridis_plasma_dark_high,
    feature_plot = TRUE,
    feature_max = 8,
    feature_col = viridis_plasma_dark_high,
    feature_raster = F,
    feature_ncol = NULL,
    version = 5,
    resolution_select = NULL,
    future.bool = F,
    ...) {
    # 0 check the settings
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(version %in% c(4, 5))
    assay_use <<- ifelse(exists("assay_use"), assay_use %||%
        DefaultAssay(object), DefaultAssay(object))
    resolution <- resolution_select %||% 0.3
    ident_group <- ident_group %||% paste0(assay_use, "_snn_res.", resolution)
    assertthat::assert_that(ident_group %in% colnames(object@meta.data))
    feature_ncol <- feature_ncol %||% round(sqrt(feature_max))
    parallel.bool <- ifelse(future.bool == T, "future.apply", "None")
    # 1 clean the markers
    if (!is.list(markers)) {
        markers <- as.list(markers)
        message("The markers should be a list object whose names are celltypes!")
    }
    scaled_markers <- lapply(markers, function(m) {
        check_markers(m, object)
    }) %>% list_clean()
    assertthat::assert_that(length(scaled_markers) >= 1)
    dot_markers <- list_shorten(scaled_markers, dot_max) %>% list_flat()
    feature_markers <- list_shorten(scaled_markers, feature_max) %>% list_flat()
    # 2 design the plot function
    draw_dot_plot_v4 <- function(object, features, cluster, colors_use = dot_col, ident_group, ...) {
        if (is.null(cluster)) {
            DotPlot_scCustom(object,
                features = features,
                group.by = ident_group,
                flip_axes = TRUE,
                dot.scale = 20,
                colors_use = colors_use, ,
                ...
            ) + theme(axis.text.x = element_text(size = 32, angle = 90),axis.text.y = element_text(size = 32))
        } else {
            Clustered_DotPlot(object,
                features = features,
                k = as.numeric(cluster),
                plot_km_elbow = FALSE, ...
            )
        }
    }
    draw_dot_plot_v5 <- function(object, features, cluster, colors_use = dot_col, ident_group, ...) {
        # because the temporary bug of scCustomize
        if (is.null(cluster)) {
            DotPlot(object,
                features = features,
                group.by = ident_group,
                # cols = colors_use,
                ...
            ) + scale_color_gradientn(colors = colors_use) +
                theme(
                    axis.text.x = element_text(size = rel(1.4), angle = 90),
                    axis.text.y = element_text(size = rel(1.4))
                )
        } else {
            stop("need to be ")
        }
    }
    draw_feature_plot_v4 <- function(object, features, colors_use = feature_col, raster = feature_raster, feature_ncol, ...) {
        FeaturePlot_scCustom(object,
            features = features,
            colors_use = colors_use,
            raster = feature_raster,
            num_columns = feature_ncol, ...
        ) & NoAxes() & theme(plot.title = element_text(size = 40))
    }
    draw_feature_plot_v5 <- function(object, features, colors_use = feature_col, raster = feature_raster, feature_ncol, ...) {
        FeaturePlot(object,
            features = features,
            # cols = colors_use,
            raster = feature_raster,
            ncol = feature_ncol,
            order = T,
            min.cutoff = 0, ...
        ) & NoAxes() & scale_color_gradientn(colors = colors_use)
    }
    Annotation_plot <- function(plot, cell_p) {
        library(patchwork)
        plot + plot_annotation(paste0(cell_p),
            theme = theme(
                plot.title = element_text(size = 50, face = "bold")
            )
        )
    }
    spacer_plot <- function(plot, n, max) {
        library(patchwork)
        for (i in (n + 1):max) {
            plot <- plot + plot_spacer()
        }
        return(plot)
    }
    draw_dot_plot <- switch(as.character(version),
        "5" = draw_dot_plot_v5,
        "4" = draw_dot_plot_v4
    )
    draw_feature_plot <- switch(as.character(version),
        "5" = draw_feature_plot_v5,
        "4" = draw_feature_plot_v4
    )
    # 3 draw the plot
    p_dim <- DimPlot(object,
        group.by = ident_group,
        label = T,
        label.size = 10,
        pt.size = 2
    ) & NoAxes() & NoLegend() & theme(plot.title = element_text(size = 40))
    if (dot_plot == T) {
        lapply_par(seq_along(dot_markers), function(m) {
            p_dot <- draw_dot_plot(
                object = object,
                cluster = dot_cluster,
                features = dot_markers[[m]],
                colors_use = dot_col,
                ident_group = ident_group
            ) %>%
                Annotation_plot(., cell_p = names(dot_markers)[m])
            # print(p_dot)
        }, parallel = parallel.bool) %>% lapply(., print)
    }
    message("plot the dot_plot of markers")
    if (feature_plot == T) {
        lapply_par(seq_along(feature_markers), function(m) {
            p_feature <- draw_feature_plot(
                object = object,
                features = feature_markers[[m]],
                colors_use = feature_col,
                feature_ncol = feature_ncol
            ) #%>%
                #spacer_plot(plot = ., n = length(feature_markers[[m]]), max = dot_max) #+
            n = length(feature_markers[[m]])
            if (n > 6){
                design <- c(area(1,6,6,11),area(1,1,5,5))
            } else if(n < 7 && n>3) {
                design <- c(area(1,6,6,11),area(1,1,5,5),area(6,1,6,5),area(5,6,6,11))
            } else {
                design <- c(area(1,6,6,11),area(1,1,2,5),area(6,1,6,5),area(3,6,6,11))
            }
            p_combine <- p_dim  + p_feature +
                plot_layout(
                    design = design
                )
            #p_feature <- p_dim | p_feature
            p_combine %>%
                Annotation_plot(., cell_p = names(feature_markers)[m])
            # print(p_feature)
        }, parallel = parallel.bool) %>% lapply(., print)
    }
    message("plot the dot_plot of markers")
}

# plot_function <- function(markers, annotation) {
#    draw_dot_plot <- function() {
#        n_markers <- length(markers)
#        n_plots <- ceiling(n_markers / max_markers)
#        for (i in 1:n_plots) {
#            start <- (i - 1) * max_markers + 1
#            end <- min(i * max_markers, n_markers)
#            if (is.numeric(cluster)) {
#                scaled_markers <- switch(assay_use,
#                    "RNA" =  check_markers(markers, object[["RNA"]]@scale.data),
#                    "SCT" =  check_markers(markers, object[["SCT"]]@scale.data)
#                )
#                Clustered_DotPlot(object,
#                   features = scaled_markers[start:end], k = cluster,
#                    plot_km_elbow = FALSE, ...
#                )
#            } else {
#                p_dot <- DotPlot_scCustom(object, features = markers[start:end], flip_axes = TRUE, colors_use = col_dot, ...)
#            }
#            if (annotation == TRUE) {
#                print(p_dot + plot_annotation(paste0(cell_type), theme = theme(
#                    plot.title = element_text(size = 18, face = "bold")
#                )))
#            } else {
#                print(p_dot)
#            }
#        }
#    }

#    draw_feature_plot <- function() {
#        for (j in markers) {
#            p_fea <- FeaturePlot_scCustom(object,
#                features = j, colors_use = col_fea,
#                ...
#            ) & NoAxes()
#            if (annotation == TRUE) {
#                print(p_fea + plot_annotation(paste0(cell_type), theme = theme(
#                    plot.title = element_text(size = 18, face = "bold")
#                )))
#            } else {
#                print(p_fea)
#            }
#        }
#    }
#    if (dot_plot == TRUE) {
#        draw_dot_plot()
#    }

#    if (feature_plot == TRUE) {
#        draw_feature_plot()
#    }
# }

# If the unfiltered_markers is a data frame, we have multiple cell types, so we'll loop through them.
# if (class(markers) == "list") {
# Load the patchwork package so we can annotate plots
#    library(patchwork)
# Get the cell type names from the column names.
#    cell_types <- names(markers)
# Loop through the cell types
#    for (cell_type in cell_types) {
# Plot the marker genes for each cell type.
#        cat("Ploting markers in ", cell_type, "...\n")
#        plot_function(check_markers(markers[, cell_type] %>%
#            na.omit() %>% unlist() %>% as.character(), object), annotation = TRUE)
#    }
# } else {
# Otherwise, we have a single cell type, so just plot the marker genes.
#    cat("Ploting markers...\n")
#    plot_function(check_markers(na.omit(markers), object), annotation = FALSE)
# }






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
