######################################
# 1.pheatmap
# 2.ggplot2
# 3.complexpheatmap
#####################################

range_reduce <- function(df, threshold = NULL) {
        mat <- as.matrix(df)
        range_mat <- range(mat)
        cat("The range of this data is: ", range_mat, "\n")
        min_abs_range_mat <- threshold %||% min(abs(range_mat))
        mat[mat < -min_abs_range_mat] <- -min_abs_range_mat
        mat[mat > min_abs_range_mat] <- min_abs_range_mat
        return(as.data.frame(mat))
}


#-----pheatmap-------#
library(pheatmap)
pdf(paste0(dir, "/output/      .pdf"), height = 7, width = 5)
#---choose your colors of annotation----#
groupcolor <- c()
names(groupcolor) <- c()
anno_color <- list(group = groupcolor)
#---------------breaks------------------#
down
up
my_breaks <- c(-Inf, seq(-down, 0, length.out = 100), seq(0, up, length.out = 100), Inf)
## 似乎会报错，不能同时出现0，需要去重一下
#-------------Core color----------------#
my_palette <- colorRampPalette(colors = c("blue", "white", "red"))(50)
#---------------------------------------#

pheatmap(plot,
        show_rownames = F, show_colnames = F,
        cluster_cols = F, cluster_rows = F,
        border = F,
        breaks = my_breaks,
        annotation_row = ,
        annotation_col = ,
        annotation_colors = anno_colors,
        annotation_legend = F, annotation_names_row = F, annotation_names_col = F,
        color = my_palette,
        treeheight_row = 0, treeheight_col = 0,
        gaps_row = c(), gaps_col = c()
)
dev.off()
