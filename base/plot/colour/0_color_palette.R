

 suppressPackageStartupMessages(library(tidyverse))
 
choose_pal <- function(n,
              source=c('viridis','all','picture','cycle'),
              name,
              picture_path=NULL){
  suppressPackageStartupMessages(library(viridis))
  suppressPackageStartupMessages(library(paletteer))
  # 1 viridis主题色
  choose_pal_viridis <- function(n=n,name=name){
    ##viridis
    #viridis：option D，为默认色带，翠绿色；
    #magma：option A，岩溶色；
    #inferno：option B，火焰色；
    #plasma： option C，血色；
    #cividis： option E;
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
  # 2 paletteer全色彩 62+8个包
  choose_pal_all <- function(){
    #palettes_d_names %>% 
    #  select(package) %>% 
    #    unique %>% 
    #    print(n=nrow(.))
    #palettes_c_names %>% 
    #  select(package) %>% 
    #    unique %>% 
    #      print(n=nrow(.))    
    tryCatch({
      paletteer_c(name, n = n) %>%  fct_inorder()
    },error=function(e){
      paletteer_d(name, n = n) %>%  fct_inorder()
    })
    
  }

  # 3 图片取色
  choose_pal_picture <- function(picture_path=picture_path){
    rPlotter::extract_colours(picture_path)
  }
 
  # 4 色轮取色
  choose_pal_cycle <- function(name,n){
    NULL
  }
  
  switch(source,'viridis'=choose_pal_viridis,
                'paletteer'=choose_pal_all,
                'picture'=choose_pal_picture,
                'cycle'=choose_pal_cycle
            )
}

 
 

##RColorBrewer
#brewer.pal.info
#pal <- RColorBrewer::brewer.pal(n, name)
# n代表选择该调色板的几个颜色出来，name代表选择的调色板的名字










# 3 ggplot叠加
##viridis
#+scale_fill_viridis(option="magma")#连续
#+scale_fill_viridis_d(option="magma")#离散
#+ scale_color_paletteer_d("basetheme::ink")
#+ scale_color_paletteer_c("ggthemes::Red-Gold")
#+ scale_color_paletteer_c(palette, direction = 1, dynamic = FALSE, ...)



# 颜色指南
color_guide <- function(path){
  if(T){
    palettes_d_names %>% mutate(package_title=package) %>%
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
  pdf(sstr_c(path,'/color_album_d.pdf'))
    for(i in seq(1, length(p.list),10)){
      if(i+9 <= length(p.list)){
        p.list[i:(i+9)] %>% patchwork::wrap_plots() %>% print
      }else {
        p.list[i:length(p.list)] %>% patchwork::wrap_plots() %>% print
      }
    }
  dev.off()
}

 }
if(T){
palettes_c_names %>% mutate(package_title=package) %>%
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
  pdf('/public/home/caojun/module_script/3_plot/color_album_c.pdf')
    for(i in seq(1, length(p.list),10)){
      if(i+9 <= length(p.list)){
         p.list[i:(i+9)] %>% patchwork::wrap_plots() %>% print
      }else {
         p.list[i:length(p.list)] %>% patchwork::wrap_plots() %>% print
      }
    }
  dev.off()
}

 

