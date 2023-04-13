# report the job
## when you report a long script
log_start <- function(...){
    start_time <<- Sys.time()
    cat("This job starts at:", format(start_time), "\n")
    script_path <- getwd()
    cat("This script is in:", script_path, "\n")
}
log_start()

log_done <- function(...){
    end_time <- Sys.time()
    time_difference <- end_time - start_time
    cat("This job ends at:", format(end_time), "\n")
    cat("Total time:", format(time_difference), "\n")
    print(sessionInfo())
}
## when you report a code block
#expr <- expression({...})
log_report <- function(expr,report=T) {
    if(report){
        start_time <- Sys.time()
        cat("This job starts at:", format(start_time), "\n") 
        result <- eval(expr)
        log_done()
        return(result)
    }

}