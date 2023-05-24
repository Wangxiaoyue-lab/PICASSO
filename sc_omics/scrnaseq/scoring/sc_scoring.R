sc_score_function <- function(...) {
    cat("# function list\n")
    cat("#-- sc_score_seurat: scoring by seurat\n")
    cat("#-- sc_score_aucell: scoring by aucell\n")
    cat("#-- sc_score_ucell: scoring by ucell\n")
    cat("#-- sc_score_maya: scoring by maya\n")
    cat("#-- sc_score_jasmine: scoring by jasmine\n")
}

# Seurat - AddModuleScore
sc_score_seurat <- function(object, genes_list) {
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    # 2 run
    object <- Seurat::AddModuleScore(object,
        features = genes_list,
        ctrl = 100,
        name = paste0(names(genes_list), "_seurat")
    )
    return(object)
}

# AUCell
sc_score_aucell <- function(object, genes_list) {
    library(AUCell)
    set.seed(1)
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    # 2 run
    matrix <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "data"
    ) %>% t()
    cells_rankings <- AUCell_buildRankings(matrix, nCores = 1, plotStats = TRUE)
    lapply(seq_along(genes_list), function(g_l) {
        cells_AUC <- AUCell_calcAUC(genes_list[[g_l]], cells_rankings)
        signature_exp <- data.frame(t(getAUC(cells_AUC)))
        object@meta.data[[paste0(names(genes_list)[g_l], "_aucell")]] <- signature_exp[, 1]
    })
}

# UCell
sc_score_ucell <- function(object, genes_list) {
    library(UCell)
    library(GSEABase)
    library(BiocParallel)
    set.seed(1)
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    # 2 run
    object <- AddModuleScore_UCell(object,
        features = genes_list,
        assay = NULL,
        slot = "data",
        BPPARAM = NULL,
        ncores = 1,
        chunk.size = 1000
    )
    return(object)
}

# singscore
sc_score_singscore <- function(object, genes_list) {
    library(singscore)
    set.seed(1)
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")
    matrix <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "data"
    ) %>% t()
    singscore.rank <- singscore::rankGenes(as.data.frame(matrix))
    singscore.Score <- lapply(seq_along(genes_list),function(g_l){
        singscore.scores_g_l <- singscore::simpleScore(singscore.rank,
            upSet = genes_list[[g_l]],
            centerScore = F
        ) %>% 
        as.data.frame() %>%
        dplyr::select(TotalScore)
    }) %>% do.call(cbind,.)
    colnames(singscore.Score) <- paste0(names(genes_list), "_singscore")
    object@meta.data <- cbind(object@meta.data, singscore.Score)
    return(object)
}


sc_score_plage <- function(...) {
    # too slow
}
sc_score_zscore <- function(...) {
    # too slow
}
sc_score_gsea <- function(...) {
    # too slow
}
sc_score_ssgsea <- function(...) {
    # too slow
}
sc_score_gsva <- function(...) {
    # too slow
}
sc_score_vision <- function(...) {
    # too slow
}


# Pagoda2
sc_score_pagoda2 <- function(){

}

# MAYA
sc_score_maya <- function(...) {
    next
}

# JASMINE
sc_score_jasmine <- function(object, genes_list) {
    # 0 design function
    # https://github.com/NNoureen/JASMINE
    stringsAsFactors <- FALSE
    library(stringr)
    library(GSA)
    set.seed(1)
    ##################################### STRUCTURE of JASMINE ####################
    ### Function1:-  Calculating Mean Ranks for signature genes across each cell
    RankCalculation <- function(x, genes) {
        subdata <- x[x != 0] ### Removing Dropouts from single cell
        DataRanksUpdated <- rank(subdata) ### Calculating ranks of each signature gene per cell
        DataRanksSigGenes <- DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)] ### Shortling rank vector for signature genes
        CumSum <- ifelse(length(DataRanksSigGenes), mean(DataRanksSigGenes, na.rm = TRUE), 0) ### Calculating Mean of ranks for signature genes
        FinalRawRank <- CumSum / length(subdata) ### Normalizing Means by total coverage
        return(FinalRawRank)
    }
    #### Function2:- Calculating enrichment of signature genes across each cell 	(using odds ratio)
    ORCalculation <- function(data, genes) {
        GE <- data[which(rownames(data) %in% genes), ] ### Subsetting data for signature genes
        NGE <- data[-which(rownames(data) %in% genes), ] ### Subsetting data for non-signature genes
        SigGenesExp <- apply(GE, 2, function(x) length(x[x != 0])) ### Calculating Number of expressed Signature Genes per cell
        NSigGenesExp <- apply(NGE, 2, function(x) length(x[x != 0])) ### Calculating Number of expressed Non-Signature Genes per cell
        SigGenesNE <- nrow(GE) - SigGenesExp ### Calculating Number of Not expressed Signature Genes per cell
        SigGenesNE <- replace(SigGenesNE, SigGenesNE == 0, 1) ### Replacing Zero's with 1
        NSigGenesExp <- replace(NSigGenesExp, NSigGenesExp == 0, 1) ### Replacing Zero's with 1
        NSigGenesNE <- nrow(data) - (NSigGenesExp + SigGenesExp) ### Calculating Number of Not expressed Non-Signature Genes per cell
        NSigGenesNE <- NSigGenesNE - SigGenesNE
        OR <- (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp) ### Calculating Enrichment (Odds Ratio)
        return(OR)
    }
    #### Function3:- Calculating enrichment of signature genes across each cell (using Likelihood ratio)
    LikelihoodCalculation <- function(data, genes) {
        GE <- data[which(rownames(data) %in% genes), ]
        NGE <- data[-which(rownames(data) %in% genes), ]
        SigGenesExp <- apply(GE, 2, function(x) length(x[x != 0]))
        NSigGenesExp <- apply(NGE, 2, function(x) length(x[x != 0]))
        SigGenesNE <- nrow(GE) - SigGenesExp
        SigGenesNE <- replace(SigGenesNE, SigGenesNE == 0, 1)
        NSigGenesExp <- replace(NSigGenesExp, NSigGenesExp == 0, 1)
        NSigGenesNE <- nrow(data) - (NSigGenesExp + SigGenesExp)
        NSigGenesNE <- NSigGenesNE - SigGenesNE
        LR1 <- SigGenesExp * (NSigGenesExp + NSigGenesNE)
        LR2 <- NSigGenesExp * (SigGenesExp + SigGenesNE)
        LR <- LR1 / LR2
        return(LR)
    }
    ###  Function 4:- Scalar [0,1] Normalization of Means and Enrichment across set of cells
    NormalizationJAS <- function(JAS_Scores) {
        JAS_Scores <- (JAS_Scores - min(JAS_Scores)) / (max(JAS_Scores) - min(JAS_Scores))
        return(JAS_Scores)
    }
    ### Function 5:- Signature Scoring via JASMINE mergining Means and Enrichment
    JASMINE <- function(data, genes, method) {
        idx <- match(genes, rownames(data))
        idx <- idx[!is.na(idx)]
        if (length(idx) > 1) {
            RM <- apply(data, 2, function(x) RankCalculation(x, genes)) ### Mean RankCalculation for single cell data matrix
            RM <- NormalizationJAS(RM) ### Normalizing Mean Ranks

            if (method == "oddsratio") {
                OR <- ORCalculation(data, genes) ### Signature Enrichment Calculation for single cell data matrix (OR)
                OR <- NormalizationJAS(OR) ### Normalizing Enrichment Scores (OR)
                JAS_Scores <- (RM + OR) / 2
            } else if (method == "likelihood") {
                LR <- LikelihoodCalculation(data, genes) ### Signature Enrichment Calculation for single cell data matrix  (LR)
                LR <- NormalizationJAS(LR) ### Normalizing Enrichment Scores (LR)
                JAS_Scores <- (RM + LR) / 2
            }
            FinalScores <- data.frame(names(RM), JAS_Scores) ### JASMINE scores
            colnames(FinalScores)[1] <- "SampleID"
            return(FinalScores)
        }
    }
    # 1 check
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(class(genes_list) == "list")

    # 2 extract the data
    data_f <- FetchData(object,
        vars = rownames(object),
        cells = colnames(object),
        layer = "counts"
    ) %>% t()

    # 3 run the program
    df_JASMINE <- lapply(seq_along(genes_list), function(g_l) {
        oddsratio <- JASMINE(data_f, genes_list[[g_l]], method = "oddsratio")
        likelihood <- JASMINE(data_f, genes_list[[g_l]], method = "likelihood")
        res_JASMINE <- cbind(oddsratio[, 2], likelihood[, 2])
        colnames(res_JASMINE) <- c(paste0(names(genes_list)[g_l], "_jasmine_or"), paste0(names(genes_list)[g_l], "_jasmine_ll"))
        return(res_JASMINE)
    }) %>%
        do.call(cbind, .) %>%
        as.data.frame()
    object@meta.data <- cbind(object@meta.data, df_JASMINE)
    return(object)
}
