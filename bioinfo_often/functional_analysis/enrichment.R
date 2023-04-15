
# do GO

profile_cluster <- function(genes = NULL, enrich_fun = "GO", species = "mm", ...) {
    library(clusterProfiler)
    Org <- switch(species,
        "mm" = "org.Mm.eg.db" %>% library(char = .),
        "hs" = "org.Hs.eg.db" %>% library(char = .)
    )
    data <- switch(enrich_fun,
        "GO" = enrichGO(gene = genes, OrgDb = names(Org), keyType = "SYMBOL", ont = "BP", ...) %>%
            .@result,
        "KEGG" = enrichKEGG(gene = genes, organism = switch(species,
            "mm" = "mmu",
            "hs" = "hsa"
        ), ...) %>%
            setReadable(OrgDb = names(Org), keyType = "ENTREZID") %>% .@result,
        "Reactome" = Reactome_function(),
        "WikiPathways" = WikiPathways_function()
    )
}


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
