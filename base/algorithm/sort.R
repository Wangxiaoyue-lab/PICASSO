sort_RRA <- function(list) {
    library(RobustRankAggreg)
    res <- aggregateRanks(list)
    return(res)
}
