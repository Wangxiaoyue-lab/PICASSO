#' Summarize pseudobulk data
#'
#' @param object An object containing data to be summarized
#' @param features A list or character vector specifying the features to summarize (default: NULL)
#' @param ... Additional arguments
#'
#' @return A data frame containing the summarized pseudobulk data
#'
#' @export
#'
summarise_pseudobulk <- function(object, features = NULL, ...) {
    UseMethod("summarise_pseudobulk", object = features)
}

summarise_pseudobulk.default <- function(object, features = NULL, ...) {
    stat_fun <- function(AE_result) {
        AE_result <- as.data.frame(AE_result)
        raw_cols <- colnames(AE_result)
        AE_result <- rownames_to_column(AE_result) %>%
            rowwise() %>%
            mutate(
                row_sum = sum(c_across(all_of(raw_cols))),
                row_mean = mean(c_across(all_of(raw_cols))),
                row_sd = sd(c_across(all_of(raw_cols))),
                row_mad = mad(c_across(all_of(raw_cols)))
            ) %>%
            ungroup() %>%
            mutate(across(all_of(raw_cols), list(norm = ~ (.x - row_mean) / row_sd), .names = "norm_{col}"))
        AE_result
    }
    result <- AverageExpression(object, features = na.omit(features), ...) %>%
        stat_fun()
    return(result)
}

summarise_pseudobulk.list <- function(object, features = NULL, ...) {
    classed_genes <- names(features)
    result <- classed_genes %>%
        map(~ {
            AE_result <- summarise_pseudobulk.default(object, features[[.x]], ...)
            AE_result$class <- .x
            AE_result
        }) %>%
        bind_rows()
    return(result)
}
