require(dplyr,quietly=T)
picasso_path <- getwd()

choose_pipeline <- function(pipeline=NULL,
                            module=NULL) {
    if(pipeline=NULL){
        list_pipeline() %>% return
    }
    if(module=NULL){
        list_module(pipeline) %>% return
    }
    load_necessary()
    lapply(pipeline,function(p){
        switch(p,
            'plot' = ifesle(module=='plot',
                    load_script(dir='visualization/plot',script='\\.R'),
                    list_module(pipeline)),
            'scrnaseq' = lapply(module,function(m){
                    switch(m,
                        'all' = load_script(dir='sc_omics/scRNA_analysis',script='\\.R'),
                        'upstream' = load_script(dir='sc_omics/scRNA_analysis/00_upstream',script='\\.R'),
                        'downstream' = load_script(dir='sc_omics/scRNA_analysis/01_modificated_seurat',script='\\.R'),
                        'pseudotime' = load_script(dir='sc_omics/scRNA_analysis/02_pseudotime',script='\\.R'),
                        'grn' = load_script(dir='sc_omics/scRNA_analysis/03_GRN',script='\\.R'),
                        'sccs' = load_script(dir='sc_omics/scRNA_analysis/04_sccs',script='\\.R'),
                        'cnv' = load_script(dir='sc_omics/scRNA_analysis/05_cnv',script='\\.R'),
                        'scoring' = load_script(dir='sc_omics/scRNA_analysis/06_scoring',script='\\.R'),
                        'chat' = load_script(dir='sc_omics/scRNA_analysis/07_chat',script='\\.R')
                        )
            })
                ,
            'bulkrnaseq' = ifesle(module=='bulkrnaseq',
                        load_script(dir='bulk_omics/bulk_rnaseq',script='\\.R'),
                        list_module(pipeline))
        )
    })
}

list_module <- function(pipeline){
    if(pipeline=='scrnaseq'){
        cat('#----scrnaseq----#\n')
        c(' --> upstream\n',
           '--> downstream\n',
           '--> pseudotime\n',
           '--> grn\n',
           '--> sccs\n',
           '--> cnv\n',
           '--> scoring\n',
           '--> chat\n',
           '--> ...\n') %>% cat
    }elseif(pipeline=='machine learning') {
        cat('#----machine learning----#\n')
        c(' --> classification\n',
           '--> regression\n',
           '--> clustering\n',
           '--> dimensionality reduction\n',
           '--> matrix decomposition\n',
           '--> ...\n') %>% cat
    }else{
        cat('the module should be just the pipeline. Or the module has not been designed.')
    }

}

list_pipeline <- function(...){
    cat('#----BASE----#\n')
    c(' --> plot\n',
       '--> ...\n',
       '--> ...\n',
       '--> ...\n',
       '--> ...\n') %>% cat 
    cat('#----sc_omics----#\n')
    c(' --> scrnaseq\n',
       '--> smartseq\n',
       '--> scatacseq\n',
       '--> spatial\n',
       '--> vdj\n') %>% cat
    cat('#----bulk_omics----#\n')
    c(' --> bulkrnaseq\n',
       '--> chipseq\n',
       '--> wes\n',
       '--> ...\n',
       '--> ...\n') %>% cat
    cat('#----stastics----#\n')
    c(' --> machine learning\n',
       '--> deep learning\n',
       '--> clinical analysis\n',
       '--> experiment analysis\n',
       '--> ...\n') %>% cat
    cat('#----other_bioinfo----#\n')
    c(' --> enrichment analysis\n',
       '--> sequence analysis\n',
       '--> annotation\n',
       '--> ...\n',
       '--> ...\n') %>% cat
}
 
# load necessary
load_necessary <- function(...){
    ## basical utils
    load_script(dir='utils',script='core_utils')
    load_script(dir='utils/parallel',script='parallel')
    load_script(dir='utils/input_your_parameter',script='parameter')
    ## color
    load_script(dir='visualization/colour',script='palette')
    ## picture
    load_script(dir='visualization/plot',script='themes')
}


# load the specified script
load_script <- function(dir,script){
    scripts <- list.files(path = paste0(picasso_path,'/',dir),
               pattern = paste(script,collapse='|'),
               recursive=T,full=T)
    lapply(scripts,source)
}




