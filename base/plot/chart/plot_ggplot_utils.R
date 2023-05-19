# 直方图x轴limits范围
calc_x_limits <- function(data) {
    p <- ggplot(data.frame(data), aes(x = data)) +
        geom_histogram(binwidth = 0.4)
    hist_data <- ggplot_build(p)$data[[1]]
    peak <- hist_data$x[which.max(hist_data$count)]
    x_max <- max(hist_data$x)
    x_min <- min(hist_data$x)
    dist <- max(abs(peak - x_max), abs(peak - x_min))
    c(peak - dist - 5, peak + dist + 5)
}


# 输入x的limits范围得到breaks
get_breaks <- function(data) {
    xmin <- min(data)
    xmax <- max(data)
    range <- ceiling(xmax) - round(xmin)
    if (range <= 10) {
        by <- 1
    } else if (range <= 40) {
        by <- 5
    } else if (range <= 60) {
        by <- 10
    } else if (range <= 200) {
        by <- 15
    } else {
        by <- 30
    }
    breaks <- seq(floor(xmin / by) * by, ceiling(xmax / by) * by, by = by)
    return(breaks)
}
