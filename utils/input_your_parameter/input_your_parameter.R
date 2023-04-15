library(optparse)
## reading parameters
option_list <- list(  
  make_option(c("-c", "--count"), type="character", 
              help="count file path"),
  make_option(c("-m", "--meta"), type="character", 
              help="meta (celltype) file path"),
  make_option(c("-o", "--output"), type="character",
              help="output dir"),
  make_option(c("-d", "--delim"), type="character", default = '@',help="the delimiter used to separate barcode and celltype in cell label"),
  make_option(c("-f", "--field"), type="integer", default = 2,help="the field of celltype in the cell label")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

count_file <- opts$count
meta_file <- opts$meta
output_dir <- opts$output
names.delim <- opts$delim
names.field <- opts$field


