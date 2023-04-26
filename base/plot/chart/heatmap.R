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
                        annotation_row = anno_col,
                        annotation_col = anno_row,
                        annotation_colors = anno_colors,
                        annotation_legend = F, annotation_names_row = F, annotation_names_col = F,
                        color = colorRampPalette(colors = palette)(100), ...
                        # treeheight_row = 0, treeheight_col = 0,
                        # gaps_row = c(), gaps_col = c()
                )
        }
}

# funkyheatmap https://github.com/funkyheatmap/funkyheatmap
