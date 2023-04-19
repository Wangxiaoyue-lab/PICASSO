require(rlang)

# library(pbapply)  #verbose for lapply
# library(parallel)
# library(doParallel)
# library(foreach)
# library(future)
# library(future.apply)
# library(snowfall)
# library(furrr)

# for Seurat
para_for_seurat <- F
if (para_for_seurat) {
    library(future)
    plan("multisession")
}

# implement a function of one variable in parallel
lapply_par <- function(x, fun,
                       parallel = c("foreach", "parallel", "snowfall", "future.apply", "None"),
                       export = NULL,
                       ncores = NULL,
                       verbose = F) {
    cl.cores <- parallel::detectCores()
    ncores <- ncores %||% cl.cores - 1
    if (parallel == "None" & verbose == TRUE) {
        # just have process bar
        return(pbapply::pblapply(X = x, FUN = fun))
    }
    if (parallel == "None") {
        # just nothing
        return(base::lapply(X = x, FUN = fun))
    }
    if (parallel == "parallel") {
        # easy for user but hard for memory
        library(parallel)
        cl <- makeCluster(getOption("cl.cores", ncores))
        clusterExport(cl, "x", "fun", export)
        par_res <- parallel::parLapply(cl, X = x, FUN = fun)
        stopCluster(cl)
        return(par_res)
    }
    if (parallel == "foreach") {
        # need doParallel thatr is based on foreach,iterators and parallel
        library(doParallel)
        library(foreach)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        outlist <- foreach::foreach(i = 1:length(x), .export = export) %dopar% {
            fun(x[[i]])
        }
        stopCluster(cl)
        names(outlist) <- names(x)
        return(outlist)
    }
    if (parallel == "snowfall") {
        library(snowfall)
        sfInit(parallel = TRUE, cpus = ncores)
        sfExport("X", "fun")
        par_res <- snowfall::sfLapply(X = x, fun)
        sfStop()
        return(par_res)
    }
    if (parallel == "future.apply") {
        library(future)
        library(future.apply)
        options(future.globals.maxSize= 1e10)
        plan("multisession", workers = ncores)
        return(future.apply::future_lapply(X = x, FUN = fun))
    }
}

# implement a nested loop function in parallel
# maybe wrong!
lapply_nest <- function(x, y, fun, parallel = c("foreach", "future.apply"), ncores = NULL) {
    cl.cores <- parallel::detectCores()
    ncores <- ncores %||% cl.cores - 1
    if (parallel == "foreach") {
        library(doParallel)
        library(foreach)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        res <- foreach(x, .combine = "c") %.%
            foreach(y) %dopar% {
                fun(x, y)
            }
    }
    if (parallel == "future.apply") {
        library(future)
        library(future.apply)
        plan("multisession", works = nocres)
        res <- future_lapply(x, function(y) {
            future_lapply(y, fun)
        })
    }
}


map_par <- function(x, fun, ncores = NULL) {
    cl.cores <- parallel::detectCores()
    ncores <- ncores %||% cl.cores - 1
    plan("multisession", works = nocres)
    return(future_map(x, fun))
}







if (F) {
    system.time({
        test <- lapply_par(1:100000, cumsum, parallel = "None")
    })
    system.time({
        test <- lapply_par(x = 1:1000000, fun = cumsum, parallel = "parallel")
    })
    system.time({
        test <- lapply_par(1:1000000, cumsum, parallel = "foreach")
    })
    system.time({
        test <- lapply_par(1:100000, cumsum, parallel = "None", verbose = T)
    })
    system.time({
        test <- lapply_par(1:1000000, cumsum, parallel = "future.apply")
    })
}
