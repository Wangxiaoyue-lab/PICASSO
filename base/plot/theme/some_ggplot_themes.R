require(ggplot2)

#空白极简
blank_all <- theme_bw()+theme(
    panel.grid=element_blank(),
    panel.background=element_rect(fill='white',color='black'),
    axis.txt=element_blank(),
    axis.ticks=element_blank(),
    panel.border=element_blank()
)

# 背景
theme_background <- theme(panel.background = element_rect(fill = 'grey80'),
          plot.background=element_rect(fill="khaki"), 
          plot.margin = unit(c(3, 2, 1, 1), "cm"))

#字号
theme_font_size <- theme(
        axis.text=element_text(size=rel(size)),
        strip.text=element_text(size=rel(size)),
        legend.text=element_text(size=rel(size)),
        plot.title=element_text(size=rel(size))
    )
 

# 坐标轴
theme_axis <- theme(axis.title.x=element_text(vjust=1,  
                                      size=20),  # X axis title
            axis.title.y=element_text(size=10,
                                     color = "blue"),  # Y axis title
            axis.text.x=element_text(size=10, 
                                     angle = 45,
                                     color = "red",
                                     vjust=.5),  # X axis text
            axis.text.y=element_text(size=10))  # Y axis text


# 图例
theme_legend <- theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
           legend.justification=c(1,0), 
           legend.position=c(0.95, 0.05),  
           legend.background = element_blank(),
           legend.key = element_blank())
