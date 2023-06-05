######################################
# 1.pheatmap
# 2.ggplot2
# 3.complexpheatmap
#####################################

#' Adjust the range of a data frame
#'
#' This function takes a data frame and an optional threshold value. It adjusts the range of the data frame so that all values are within the range of [-threshold, threshold]. If no threshold is provided, the function uses the minimum absolute value of the range of the data frame as the threshold.
#'
#' @param df A data frame.
#' @param threshold An optional numeric value specifying the threshold for adjusting the range of the data frame.
#' @return A data frame with adjusted range.
#' @examples
#' df <- data.frame(x = c(-10, 0, 10), y = c(-5, 0, 5))
#' plot_range_adjust(df)
#' plot_range_adjust(df, threshold = 3)
plot_range_adjust <- function(df, threshold = NULL) {
        mat <- as.matrix(df)
        range_mat <- range(mat)
        cat("The range of this data is: ", range_mat, "\n")
        min_abs_range_mat <- threshold %||% min(abs(range_mat))
        mat[mat < -min_abs_range_mat] <- -min_abs_range_mat
        mat[mat > min_abs_range_mat] <- min_abs_range_mat
        return(as.data.frame(mat))
}


plot_gap <- function(vector) {
        which(!duplicated(vector))[-1] - 1
}


plot_rearrange <- function(df, group_df, seed = 1) {
        set.seed(seed)
        group <- group_df %>%
                dplyr::mutate(names = rownames(.)) %>%
                dplyr::arrange(across(-names)) %>%
                group_by(across(-names)) %>%
                nest()
        order <- lapply(group$data, function(d) {
                tryCatch(
                        {
                                df[, d %>% pull(names)] %>%
                                        t() %>%
                                        dist() %>%
                                        hclust() %>%
                                        .$order
                        },
                        error = function(e) {
                                seq_along(d %>% pull(names))
                        }
                )
        })
        arrange <- mapply(function(x, y) {
                x %<>% pull(names)
                x[y]
        }, group$data, order, SIMPLIFY = T) %>% unlist()
        df <- df[, match(arrange, colnames(df))]
        return(df)
}


#' Reorder hierarchical clustering by the specified method
#'
#' This function reorders hierarchical clustering of a matrix based on the specified method. It supports several methods, such as "raw", "manual", "row_select", "col_names", and "eigenvalue".
#'
#' @param mat a matrix or data frame for clustering
#' @param method character string specifying the reordering method, one of "raw", "manual", "row_select", "col_names", "eigenvalue"
#' @param manual optional vector specifying the manual order for rows, only used if method is "manual"
#' @param row_select optional character string specifying a column in mat for row ordering, only used if method is "row_select"
#' @return an object of class `hclust` with reordered clustering
#' @import dendextend
#' @import assertthat
#' @export
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10)
#' ordered_hclust <- hclust_order(mat, method = "raw")
hclust_order <- function(mat, method, manual = NULL, row_select = NULL) {
        library(dendextend)
        assertthat::assert_that(method %in% c("raw", "manual", "row_select", "col_names", "eigenvalue"))
        mat_t <- mat %>%
                t() %>%
                as.data.frame()
        hclust_mat <- mat_t %>%
                dist() %>%
                hclust()
        if (method == "raw") {
                return(hclust_mat)
        } else if (method == "manual") {
                hclust_manual <- reorder(as.dendrogram(hclust_mat),
                        wts = order(match(manual, rownames(mat_t)))
                ) %>%
                        as.hclust()
                return(hclust_manual)
        } else if (method == "row_select") {
                hclust_row_select <- reorder(as.dendrogram(hclust_mat),
                        wts = mat_t[, row_select]
                ) %>%
                        as.hclust()
                return(hclust_row_select)
        } else if (method == "col_names") {
                hclust_col_names <- hclust_mat %>%
                        as.dendrogram() %>%
                        sort() %>%
                        as.hclust()
                return(hclust_col_names)
        } else if (method == "eigenvalue") {
                hclust_eigenvalue <- reorder(as.dendrogram(hclust_mat),
                        wts = svd(mat)$v[, 1]
                ) %>%
                        as.hclust()
                return(hclust_eigenvalue)
        }
}



#' Plot a heatmap
#'
#' This function takes a data matrix and optional parameters to plot a heatmap using the pheatmap package. The function allows for customization of the color palette, row and column annotations, clustering options, and more.
#'
#' @param data A numeric matrix or data frame containing the data to be plotted.
#' @param palette An optional character vector specifying the color palette to use for the heatmap.
#' @param anno_col An optional data frame specifying column annotations.
#' @param anno_row An optional data frame specifying row annotations.
#' @param anno_colours An optional list specifying colors for the annotations.
#' @param cluster A character string specifying whether to cluster rows, columns, both or none.
#' @param name A character string specifying whether to show row names, column names, both or none.
#' @return A heatmap plot.
#' @examples
#' data <- matrix(rnorm(100), ncol = 10)
#' plot_heatmap(data)
plot_heatmap <- function(data,
                         palette = NULL,
                         anno_col = NULL,
                         anno_row = NULL,
                         anno_colours = NULL,
                         cluster = "both", # c('both','row','col','none')
                         name = "both", # c('both','row','col','none')
                         ...) {
        library(pheatmap)
        cluster_bool <- switch(cluster,
                "both" = c(T, T),
                "row" = c(T, F),
                "col" = c(F, T),
                "none" = c(F, F)
        )
        name_bool <- switch(name,
                "both" = c(T, T),
                "row" = c(T, F),
                "col" = c(F, T),
                "none" = c(F, F)
        )
        palette <- palette %||% choose_pal(n = 7, source = "all", name = "viridis::turbo")
        sample_color <- c(
                choose_pal(n = 10, source = "all", name = "IslamicArt::shiraz"),
                choose_pal(n = 10, source = "all", name = "khroma::soil"),
                choose_pal(n = 10, source = "all", name = "IslamicArt::shiraz"),
                choose_pal(n = 10, source = "all", name = "ggthemes::stata_s2color"),
                choose_pal(n = 10, source = "all", name = "palettetown::tentacruel")
        )
        if (!is.null(anno_colours)) {
                if (class(anno_colours) != "list") {
                        stop("The anno_colour should be a list object")
                }
        }
        if (!is.null(anno_col)) {
                if (class(anno_col) != "data.frame") {
                        stop("The anno_col should be a data.frame object")
                } else {
                        anno_col_colour <- anno_col %>%
                                as.list() %>%
                                lapply(., function(g) {
                                        length(g) %>% sample(x = sample_color, size = .)
                                })
                }
        }
        if (!is.null(anno_row)) {
                if (class(anno_row) != "data.frame") {
                        stop("The anno_row should be a data.frame object")
                } else {
                        anno_row_colour <- anno_row %>%
                                as.list() %>%
                                lapply(., function(g) {
                                        length(g) %>% sample(x = sample_color, size = .)
                                })
                }
                anno_colours <- anno_colours %||% c(anno_col_colour, anno_row_colour)
                pheatmap(data,
                        show_rownames = name_bool[1],
                        show_colnames = name_bool[2],
                        cluster_rows = cluster_bool[1],
                        cluster_cols = cluster_bool[2],
                        border = F, scale = "none",
                        # breaks = my_breaks,
                        annotation_col = anno_col,
                        annotation_row = anno_row,
                        annotation_colors = anno_colors,
                        annotation_legend = F, annotation_names_row = F, annotation_names_col = F,
                        color = colorRampPalette(colors = palette)(100), ...
                        # treeheight_row = 0, treeheight_col = 0,
                        # gaps_row = c(), gaps_col = c()
                )
        }
}




#' Add flag to a pheatmap object
#'
#' This function adds flag to a pheatmap object for the given kept labels. The repel.degree parameter controls how much space is allocated for repelling labels.
#' @param pheatmap pheatmap object to add flag to
#' @param kept.labels list of labels to add flag for
#' @param repel.degree number within [0, 1], controls label repulsion. 0: spread out labels over existing range of kept labels. 1: spread out labels over the full y-axis
#' @return updated heatmap object with added flags
#' @import ggplot2
#' @import grid
#' @export
#' @examples
#' plot_add_flag(pheatmap, c("label1", "label2"), 0.5)
plot_add_flag <- function(pheatmap,
                          kept.labels,
                          repel.degree) {
        library(ggplot2)
        library(grid)

        heatmap <- pheatmap$gtable
        new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]
        new.label$label <- ifelse(new.label$label %in% kept.labels,
                new.label$label, ""
        )
        repelled.y <- function(d, d.select, k = repel.degree) {
                strip.npc <- function(dd) {
                        if (!"unit.arithmetic" %in% class(dd)) {
                                return(as.numeric(dd))
                        }
                        d1 <- strip.npc(dd$arg1)
                        d2 <- strip.npc(dd$arg2)
                        fn <- dd$fname
                        return(lazyeval::lazy_eval(paste(d1, fn, d2)))
                }
                full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
                selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
                return(unit(
                        seq(
                                from = max(selected.range) + k * (max(full.range) - max(selected.range)),
                                to = min(selected.range) - k * (min(selected.range) - min(full.range)),
                                length.out = sum(d.select)
                        ),
                        "npc"
                ))
        }
        new.y.positions <- repelled.y(new.label$y,
                d.select = new.label$label != ""
        )
        new.flag <- segmentsGrob(
                x0 = new.label$x,
                x1 = new.label$x + unit(0.15, "npc"),
                y0 = new.label$y[new.label$label != ""],
                y1 = new.y.positions
        )
        new.label$x <- new.label$x + unit(0.2, "npc")
        new.label$y[new.label$label != ""] <- new.y.positions
        heatmap <- gtable::gtable_add_grob(
                x = heatmap,
                grobs = new.flag,
                t = 4,
                l = 4
        )
        heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
        grid.newpage()
        grid.draw(heatmap)
        invisible(heatmap)
}





# 长数据热图
# stemangiola/tidyHeatmap
# https://cran.r-project.org/web/packages/tidyHeatmap/vignettes/introduction.html


# library(tidyHeatmap)
# mtcars_tidy_groupings |>
#    group_by(vs, property_group) |>  #先是行分组，再是列分组
#    mutate(`Car name` = fct_reorder(`Car name`, `Car name`, .desc = TRUE))|>    #因子排序
# heatmap(
#       .row,
#       .column,
#       .value,
#       transform = NULL,  是否进行转化，如log1p
#       scale = "none",  还有row、col、both
#       cluster_rows = FALSE,cluster_cols = FALSE,
#       column_dend_height = unit(0.2, "cm"), row_dend_width = unit(0.2, "cm"),
#       palette_value = c("#440154FF", "#21908CFF", "#fefada"),
#       show_heatmap_legend = FALSE,
#       row_names_gp = gpar(fontsize = 7),
#       column_names_gp = gpar(fontsize = 7),
#       column_title_gp = gpar(fontsize = 7),
#       row_title_gp = gpar(fontsize = 7)
#       palette_grouping = list(
# For first grouping (vs)
#            c("#66C2A5", "#FC8D62"),
# For second grouping (property_group)
#            c("#b58b4c", "#74a6aa")),
#       row_km = 2,        #kmeans聚类划分
#       column_km = 2       #kmeans聚类划分
#       )
# split_rows(2)  |>   #根据层次聚类划分
# split_columns(2)|>

# 增加注释条
# |>   add_tile(hp, show_legend = FALSE,
#               size = unit(0.3, "cm"),
#               annotation_name_gp= gpar(fontsize = 8),
#               palette = circlize::colorRamp2(c(0, 100, 200, 300), viridis::magma(4))   #选择刻度与颜色
#               palette = circlize::colorRamp2( seq(-2, 2, length.out = 11), RColorBrewer::brewer.pal(11, "RdBu") )
# add_point(activation) |>
# add_bar(size) |>
# add_line(age)|>

# 增加标记层
# layer_point(pval < 0.05)|>

# 拼图
# library(patchwork)
# wrap_heatmap(p_heatmap) + p_ggplot


# benchmark热图
# funkyheatmap https://github.com/funkyheatmap/funkyheatmap
