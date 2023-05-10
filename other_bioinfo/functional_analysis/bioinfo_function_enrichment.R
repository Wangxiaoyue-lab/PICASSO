#' Perform over-representation analysis with optional simplification
#'
#' This function takes a list of genes and performs over-representation analysis using either the Gene Ontology (GO) or KEGG pathway databases. The function uses the clusterProfiler package to perform the analysis. The results can be optionally simplified using the simplify function from the same package.
#'
#' @param genes A character vector of gene symbols.
#' @param enrich_fun A character string specifying the type of enrichment analysis to perform ("GO" or "KEGG").
#' @param species A character string specifying the species for the analysis ("mm" for mouse or "hs" for human).
#' @param simplify_cutoff An optional numeric value specifying the cutoff for semantic similarity when simplifying the results.
#' @return A data frame containing the results of the enrichment analysis.
#' @examples
#' genes <- c("BRCA1", "BRCA2", "TP53")
#' enrich_ora(genes = genes, enrich_fun = "GO", species = "hs", simplify_cutoff = 0.7)
enrich_ora <- function(genes = NULL, enrich_fun = "GO", species = "mm", simplify_cutoff = NULL, ...) {
    library(clusterProfiler)
    Org <- switch(species,
        "mm" = "org.Mm.eg.db" %>% library(char = .),
        "hs" = "org.Hs.eg.db" %>% library(char = .)
    )
    simplify_ <- function(res, simplify_cutoff = simplify_cutoff) {
        if (!is.null(simplify_cutoff)) {
            res <- clusterProfiler::simplify(res, cutoff = as.numeric(simplify_cutoff))
        }
        return(res)
    }
    data <- switch(enrich_fun,
        "GO" = enrichGO(gene = genes, OrgDb = names(Org), keyType = "SYMBOL", ont = "BP", ...) %>% simplify_() %>%
            as.data.frame(),
        "KEGG" = enrichKEGG(gene = genes, organism = switch(species,
            "mm" = "mmu",
            "hs" = "hsa"
        ), ...) %>%
            setReadable(OrgDb = names(Org), keyType = "ENTREZID") %>% as.data.frame()
    )
}

enrich_gsea <- function() {}

enrich_gsva <- function() {}

enrich_kobas <- function() {}
# http://kobas.cbi.pku.edu.cn/download/

enrich_go_module <- function() {}


#' Plot GO enrichment results
#'
#' This function takes a data frame of GO enrichment results and creates a bar plot of the top n most significant terms. The function allows for grouping and coloring of the bars based on user-specified columns in the data frame.
#'
#' @param df A data frame containing GO enrichment results.
#' @param group An optional character string specifying the column to use for grouping the results. For example, different cell types.
#' @param class An optional character string specifying the column to use for coloring the bars. Must be "up" and/or "done".
#' @param n An integer specifying the number of top terms to plot by "p.adjust".
#' @return A bar plot of the top n most significant GO terms.
#' @examples
#' data <- enrich_ora(gene = c("BRCA1", "BRCA2", "TP53"), enrich_fun = "GO", species = "mm", simplify_cutoff = NULL)
#' GO_plot(data)
GO_plot <- function(df, group = NULL, class = NULL, n = 10) {
    # create plot function
    plot_func <- function(sub_df) {
        if (!is.null(class)) {
            sub_df <- sub_df %>%
                group_by(.[[class]]) %>%
                top_n(n, -p.adjust) %>%
                ungroup()
        } else {
            sub_df <- sub_df %>%
                top_n(n, -p.adjust)
        }
        sub_df <- sub_df %>%
            distinct(ID, .keep_all = TRUE) %>%
            mutate(log_p_adjust = ifelse(!is.null(class) & .data[[class]] == "down",
                log10(p.adjust), -log10(p.adjust)
            ))
        p <- ggplot(sub_df, aes(
            x = log_p_adjust,
            y = reorder(Description, log_p_adjust),
        )) +
            geom_col() +
            theme_bw() +
            scale_y_discrete(labels = function(x) {
                str_wrap(x, width = quantile(nchar(sub_df$Description), 0.75))
            })
        if (!is.null(class)) {
            p <- p + geom_col(aes(fill = .data[[class]])) +
                scale_fill_manual(values = c("red", "blue"), breaks = c("up", "down"))
        } else {
            p <- p + geom_col()
            +scale_fill_gradient(high = "red", low = "blue")
        }
        if (!is.null(group)) {
            p <- p + ggtitle(unique(sub_df[[group]]))
        }
        print(p + ylab("Description"))
    }

    # apply plot function
    if (is.null(group)) {
        plot_func(df)
    } else {
        df %>%
            group_by(df[[group]]) %>%
            group_map(~ plot_func(.x))
    }
}
