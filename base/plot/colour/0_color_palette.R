
 suppressPackageStartupMessages(library(viridis))
#library(RColorBrewer)
 suppressPackageStartupMessages(library(paletteer))
#library(rPlotter)
 suppressPackageStartupMessages(library(tidyverse))
# 1 主题色提供

pal_palette <- function(n,source,name)


## paletteer 62+8个包
d_palettes <- palettes_d_names
d_palettes %>% select(package) %>% unique %>% print(n=62)
c_palettes <- palettes_c_names 
c_palettes %>% select(package) %>% unique %>% print(n=8)

pal <- paletteer_c("scico::berlin", n = 10) %>%  fct_inorder()
paletteer_d("RColorBrewer::Paired")
paletteer_dynamic("cartography::green.pal", 5)


##RColorBrewer
brewer.pal.info
pal <- RColorBrewer::brewer.pal(n, name)
# n代表选择该调色板的几个颜色出来，name代表选择的调色板的名字

##viridis
#viridis：option D，为默认色带，翠绿色；
#magma：option A，岩溶色；
#inferno：option B，火焰色；
#plasma： option C，血色；
#cividis： option E;
pal_viridis <- function(n,name){
    switch(name,viridis::viridis = viridis(n),
                viridis::magma = magma(n),
                viridis::inferno = inferno(n),
                viridis::plasma = plasma(n),
                viridis::cividis = cividis(n),
                viridis::rockect = rockect(n),
                viridis::mako = mako(n),
                viridis::turbo = turbo(n)
                )
}






# 2 图片取色
pal <- rPlotter::extract_colours(picture_path)


# 3 ggplot叠加
##viridis
#+scale_fill_viridis(option="magma")#连续
#+scale_fill_viridis_d(option="magma")#离散
+ scale_color_paletteer_d("basetheme::ink")
+ scale_color_paletteer_c("ggthemes::Red-Gold")
+ scale_color_paletteer_c(palette, direction = 1, dynamic = FALSE, ...)



# 菜谱
if(F){
all_palettes <- palettes_d_names
all_palettes %>% mutate(package_title=package) %>%
  group_by(package) %>%
  nest() %>%
  mutate(
    plots = map(data, ~{
      .x %>%  #是data
        mutate(
          palette_res = map2(palette,package_title, ~{
            tryCatch(
              {
                unlist(paletteer_d(palette=paste0(.y,"::",.x)))
              },
              error = function(e) {
                NA
              }
            )
          })
        ) %>%
        filter(!map_lgl(palette, is.na)) %>%
        mutate(
          plot = pmap(list(package_title, palette, palette_res), ~{
            ggplot() +
              geom_tile(aes(x = seq_along(..3), y = 0, fill = ..3)) +
              scale_fill_identity() +
              labs(title = paste(..1, ..2)) +
              theme_void() +
              theme(legend.position="none", plot.title=element_text(size=8))
          })
        ) %>%
        pull(plot)
    })
  ) %>%
  pull(plots) %>%
  flatten() %>% as.list() -> p.list
}
pdf('/public/home/caojun/module_script/3_plot/color_album_d.pdf')
 #
  for(i in seq(1, length(p.list),10)){
    if(i+9 <= length(p.list)){
        p.list[i:(i+9)] %>% patchwork::wrap_plots() %>% print

    }else {
       p.list[i:length(p.list)] %>% patchwork::wrap_plots() %>% print
    }

  }
dev.off()


if(F){
all_palettes <- palettes_c_names
all_palettes %>% mutate(package_title=package) %>%
  group_by(package) %>%
  nest() %>%
  mutate(
    plots = map(data, ~{
      .x %>%  #是data
        mutate(
          palette_res = map2(palette,package_title, ~{
            tryCatch(
              {
                unlist(paletteer_c(palette=paste0(.y,"::",.x),n=10))
              },
              error = function(e) {
                NA
              }
            )
          })
        ) %>%
        filter(!map_lgl(palette, is.na)) %>%
        mutate(
          plot = pmap(list(package_title, palette, palette_res), ~{
            ggplot() +
              geom_tile(aes(x = seq_along(..3), y = 0, fill = ..3)) +
              scale_fill_identity() +
              labs(title = paste(..1, ..2)) +
              theme_void() +
              theme(legend.position="none", plot.title=element_text(size=8))
          })
        ) %>%
        pull(plot)
    })
  ) %>%
  pull(plots) %>%
  flatten() %>% as.list() -> p.list
}
pdf('/public/home/caojun/module_script/3_plot/color_album_c.pdf')
  for(i in seq(1, length(p.list),10)){
    if(i+9 <= length(p.list)){
        p.list[i:(i+9)] %>% patchwork::wrap_plots() %>% print

    }else {
       p.list[i:length(p.list)] %>% patchwork::wrap_plots() %>% print
    }
  }
dev.off()
 

