sc_dc_nmf <- function(object, nfeatures, k, k_downsample=TRUE,matrix_downsample=TRUE,n_top=NULL){
    #V=W*H
    # check
    assertthat::assert_that(class(object) == "Seurat")
    matrix2nmf <- object %>% NormalizeData() %>%
        FindVariableFeatures(nfeatures = nfeatures) %>%
        ScaleData(do.center = F) %>% 
        FetchData(.,
            vars = VariableFeatures(.),
            cells = colnames(.),
            layer = "scale.data"
        ) %>%
        t
    n_top <- n_top %||% 200
    if(length(k)==1){
        res_best <- NMF::nmf(matrix2nmf, k, seed = 1)
        k.best  <- k
        p_k <- NULL
    } else {
       if(k_downsample==TRUE){
            select_numbers <- function(numbers) {
                n <- length(numbers)
                if (n <= 5) {
                    result <- numbers
                } else if (n <= 20) {
                    result <- round(seq(from = numbers[1], to = numbers[n], length.out = 5))
                } else if (n <= 50) {
                    result <- round(seq(from = numbers[1], to = numbers[n], length.out = 10))
                } else if (n <= 100) {
                    result <- round(seq(from = numbers[1], to = numbers[n], length.out = 20))
                } else {
                    result <- round(seq(from = numbers[1], to = numbers[n], length.out = 30))
                }

                return(result)
            }
            k_set <- select_numbers(k)
       }else {
            k_set <- k
       }
       if(matrix_downsample==TRUE){
            matrix_nmf <- matrix2nmf[, sample(1:ncol(matrix2nmf), min(max(k) * 200), ncol(matrix2nmf) / 5)]
       }else{
            matrix_nmf <- matrix2nmf
       }
        res_pre <- NMF::nmf(matrix_nmf, k_set, seed = 1)
        coph <- res_pre$measures$cophenetic
        coph_diff <- NULL
        for (i in 2:length(coph))
        {
            coph_diff <- c(coph_diff, coph[i - 1] - coph[i])
        }
        k.best <- k_set[which.max(coph_diff) + 1]
        p_k <- plot(k_set, coph, type = "b", col = "purple")
        res_best <- NMF::nmf(matrix2nmf, k.best, seed = 1)
    }
    # 计算重构误差
    ## 计算 Frobenius 范数
    #frobenius_norm <- function(V, W, H) {
    #    return(sqrt(sum((V - W %*% H)^2)))
    #}
    # 每因子的重要特征
    top_features <- extractFeatures(res_best, as.integer(n_top)) %>%
        lapply(., function(x) {
            rownames(res_best)[x]
        }) %>%
        do.call("rbind", .)
    # 添加到object中
    object <- RunPCA(object, reduction.name = "nmf", verbose = F)
    object@reductions$nmf@cell.embeddings <- t(coef(res_best))
    object@reductions$nmf@feature.loadings <- basis(res_best)
    colnames_nmf <- paste0("NMF_cluster_", k.best)
    object@meta.data %<>%
    mutate(!!sym(colnames_nmf) = apply(NMF::coefficients(res_best), 2, which.max))
    object_list <- list(object = object, 
        top_features = top_features, 
        plot_cophenetic = p_k)
    message("#output\n1.seurat object\n2.top features\n3.plot of cophenetic")
    return(object_list)
}




sc_dc_cnmf <- function(){

}