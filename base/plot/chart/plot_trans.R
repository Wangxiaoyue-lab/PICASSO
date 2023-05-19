# from other plot object(e.g. pheatmap) to ggplot2
if (F) {
    require(ggplotify)
    gg_p <- as.ggplot(p)
}


# from pictures to ggplot2
if (F) {
    library(magick)
    require(ggplotify)
    img <- image_read(path)
    p <- as.ggplot(img)
}

# from ggplot2 to data
if (F) {
    ggplot_build(p)
}

# from ggplot2 to code
if (F) {
    require(ggreverse)
    plot_code <- ggreverse::convert_to_code(p)
    # 重新执行
    eval(parse(text = plot_code))
}
