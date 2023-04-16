# summarise
## summarise the expression of pseudobulk group
summarise_pseudobulk <- function(object, features = NULL, ...) {
    # @ features: should be `list`
    stat_fun <- function(AE_result) {
        AE_result <- as.data.frame(AE_result)
        raw_cols <- colnames(AE_result)
        AE_result <- rowwise(AE_result) %>%
            mutate(
                row_sum = sum(c_across(all_of(raw_cols))),
                row_mean = mean(c_across(all_of(raw_cols))),
                row_sd = sd(c_across(all_of(raw_cols))),
                row_mad = mad(c_across(all_of(raw_cols)))
            ) %>%
            ungroup() %>%
            mutate(
                across(all_of(raw_cols), list(norm = ~ (.x - row_mean) / row_sd), .names = "norm_{col}")
            )
    }

    if (class(features) == "list") {
        classed_genes <- names(features)
        result <- classed_genes %>%
            map(~ {
                AE_result <- AverageExpression(object, features = as.character(na.omit(features[[.x]]), ...)) %>% stat_fun()
                AE_result$class <- .x
                AE_result
            }) %>%
            bind_rows()
    } else {
        result <- AverageExpression(object, features = na.omit(features), ...) %>% stat_fun()
    }
    return(result)
}