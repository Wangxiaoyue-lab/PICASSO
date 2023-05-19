plot_pval <- function(plot, pairs = list(c(1, 2)), test = "wilcox.test") {
    library(ggpval)
    assertthat::assert_that(any(class(plot) == "ggplot"))
    assertthat::assert_that(test %in% c("wilcox.test", "t.test"))
    add_pval(plot, pairs = pairs, test = test)
}
