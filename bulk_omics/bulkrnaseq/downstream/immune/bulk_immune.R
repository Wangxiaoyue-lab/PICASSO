library(immunedeconv)
# devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)

# https://omnideconv.org/immunedeconv/
# https://omnideconv.org/immunedeconv/articles/immunedeconv.html#special-case-cibersort


# immunedeconv
## quantiseq  10
## timer  6
## cibersort  22
## cibersort_abs
## mcp_counter    8
## xcell  64
## epic   6
## abis
## consensus_tme
## estimate


immunedeconv::deconvolute(gene_expression = gene_expression, methods = "")

# indications
immunedeconv::timer_available_cancers

rna_immune <- function(data, methods, indication) {
    data_cibersort <-
        data_epic <- immunedeconv::deconvolute(data, "epic")
    data_xcell <- immunedeconv::deconvolute(data, "xcell")
    data_estimate <- immunedeconv::deconvolute(data, "estimate")
    data_quantiseq <- immunedeconv::deconvolute(data, "quantiseq")
    data_abis <- immunedeconv::deconvolute(data, "abis")
    data_mcp_counter <- immunedeconv::deconvolute(data, "mcp_counter")
    data_consensus_tme <- immunedeconv::deconvolute(data, "consensus_tme",
        indications = rep(indication, ncol(data))
    )
    data_timer <- immunedeconv::deconvolute(data, "timer",
        indications = rep(indication, ncol(data))
    )
}


# IOBR
## ips
## svr
## lsei

# https://iobr.github.io/IOBR/IOBR-VIGNETTE.html
# https://github.com/IOBR/IOBR
tme_deconvolution_methods
#>         MCPcounter               EPIC              xCell          CIBERSORT
#>       "mcpcounter"             "epic"            "xcell"        "cibersort"
#> CIBERSORT Absolute                IPS           ESTIMATE                SVR
#>    "cibersort_abs"              "ips"         "estimate"              "svr"
#>               lsei              TIMER          quanTIseq
#>             "lsei"            "timer"        "quantiseq"
# Return available parameter options of TME deconvolution.
