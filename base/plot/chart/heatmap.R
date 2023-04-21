######################################
# 1.pheatmap
# 2.ggplot2
# 3.complexpheatmap
#####################################

plot_range_adjust <- function(df, threshold = NULL) {
        mat <- as.matrix(df)
        range_mat <- range(mat)
        cat("The range of this data is: ", range_mat, "\n")
        min_abs_range_mat <- threshold %||% min(abs(range_mat))
        mat[mat < -min_abs_range_mat] <- -min_abs_range_mat
        mat[mat > min_abs_range_mat] <- min_abs_range_mat
        return(as.data.frame(mat))
}


#-----pheatmap-------#
plot_heatmap <- function(data,
                        palette=NULL,
                        anno_col=NULL,
                        anno_row=NULL,
                        anno_colours=NULL,
                        cluster='both',#c('both','row','col','none')
                        name='both',#c('both','row','col','none')
                        ...
                        ){
        library(pheatmap)
        cluster_bool <- switch(cluster,
                                'both'=c(T,T),
                                'row'=c(T,F),
                                'col'=c(F,T),
                                'none'=c(F,F))
        name_bool <- switch(name,
                                'both'=c(T,T),
                                'row'=c(T,F),
                                'col'=c(F,T),
                                'none'=c(F,F))
        palette <- palette %||% choose_pal(n=7,source='all',name='viridis::turbo')
        sample_color <- c(choose_pal(n=10,source='all',name='IslamicArt::shiraz'),
                          choose_pal(n=10,source='all',name='khroma::soil'),
                          choose_pal(n=10,source='all',name='IslamicArt::shiraz'),
                          choose_pal(n=10,source='all',name='ggthemes::stata_s2color'),
                          choose_pal(n=10,source='all',name='palettetown::tentacruel'))
        if(!is.null(anno_colours)){
                if(class(anno_colours)!='list'){
                        stop('The anno_colour should be a list object')
                }
        }
        if(!is.null(anno_col)){
                if(class(anno_col)!='data.frame'){
                        stop('The anno_col should be a data.frame object')
                }else{
                        anno_col_colour <- anno_col %>% as.list %>% lapply(.,function(g){
                                length(g) %>% sample(x=sample_color,size=.)
                        })
                }
        }
        if(!is.null(anno_row)){
                if(class(anno_row)!='data.frame'){
                        stop('The anno_row should be a data.frame object')
                }else{
                        anno_row_colour <- anno_row %>% as.list %>% lapply(.,function(g){
                                length(g) %>% sample(x=sample_color,size=.)
                        })
        }
        anno_colours <- anno_colours %||% c(anno_col_colour,anno_row_colour)
        if
        pheatmap(data,
                show_rownames = name_bool[1], 
                show_colnames = name_bool[2],
                cluster_rows = cluster_bool[1],
                cluster_cols = cluster_bool[2], 
                border = F,scale='none',
                #breaks = my_breaks,
                annotation_row = anno_col,
                annotation_col = anno_row,
                annotation_colors = anno_colors,
                annotation_legend = F, annotation_names_row = F, annotation_names_col = F,
                color = colorRampPalette(colors = palette )(100),...
                #treeheight_row = 0, treeheight_col = 0,
                #gaps_row = c(), gaps_col = c()
                )
        }
                        }

 