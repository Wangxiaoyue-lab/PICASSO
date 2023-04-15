# packages and scripts



source("../utils/core_utils.R")
source("../utils/load_ref.R")
getwd() %>%
    paste0(., "/01_modificated_seurat") %>%
    list.files(., pattern = "\\.R$", full.names = TRUE) %>%
    lapply(., source)

# set important global variables

# pal


