t_start <- Sys.time()

suppressPackageStartupMessages(
    suppressWarnings( {
        library(foreach)
        library(doParallel)
    }))
    


cat("Preparing parallel simulations... \n")

script_dir <- file.path(getwd(), "R")
data_dir <- file.path(getwd(), "data", "sim.generated")
log_dir <- file.path(getwd(), "R", "rlogs")
if(!dir.exists(log_dir)) dir.create(log_dir)

p_mtds <- c("p.seurat.R", "p.spark.R", "p.sparkx.R", "p.nnsvg.R", "p.giotto.R", "p.senna.R")
r_mtds <- c("r.seurat.R", "r.spark.R", "r.sparkx.R", "r.nnsvg.R", "r.giotto.R", "r.senna.R")
i_mtds <- c("i.seurat.R", "i.spark.R", "i.sparkx.R", "i.nnsvg.R", "i.giotto.R", "i.senna.R")
si <- seq(from = 0.1, to = 1, length.out = 30)





n_cores <- 10
cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

cat("\n=== Starting progression benchmark simulations... ===\n")

res_p <- foreach(siv = si, .combine = c) %dopar% {
    log_msg <- paste0(">>> [Worker] Processing SI: ", round(siv, 2), "\n")
    cat(log_msg)
    
    for(mtd in p_mtds) {
        cmd <- paste("taskset -c 0-11 Rscript", file.path(script_dir, mtd), siv)
        
        log_file <- file.path(log_dir, paste0(mtd, ".p.SI", round(siv, 2), ".log"))
        cmdq <- paste(cmd, ">", log_file, "2>&1")

        exit_code <- system(cmdq)
        if(exit_code != 0) {
            err_msg <- paste0("!!! [Worker] Error in method ", mtd, " at SI ", round(siv, 2), "\n")
            cat(err_msg)
        }
    }

    return(paste0("SI ", siv, " Completed."))
}

stopCluster(cl)
gc()
cat("\n Progression Done. \n")





cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

cat("\n=== Starting regionation benchmark simulations... ===\n")

res_r <- foreach(siv = si, .combine = c) %dopar% {
    log_msg <- paste0(">>> [Worker] Processing SI: ", round(siv, 2), "\n")
    cat(log_msg)
    
    for(mtd in r_mtds) {
        cmd <- paste("taskset -c 0-11 Rscript", file.path(script_dir, mtd), siv)
        
        log_file <- file.path(log_dir, paste0(mtd, ".r.SI", round(siv, 2), ".log"))
        cmdq <- paste(cmd, ">", log_file, "2>&1")

       exit_code <- system(cmdq)
       if(exit_code != 0) {
           err_msg <- paste0("!!! [Worker] Error in method ", mtd, " at SI ", round(siv, 2), "\n")
           cat(err_msg)
       }
   }

   return(paste0("SI ", siv, " Completed."))
}

stopCluster(cl)
gc()
cat("\n Regionation Done. \n")





cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

cat("\n=== Starting islet benchmark simulations... ===\n")

res_i <- foreach(siv = si, .combine = c) %dopar% {
    log_msg <- paste0(">>> [Worker] Processing SI: ", round(siv, 2), "\n")
    cat(log_msg)
    
    for(mtd in i_mtds) {
        cmd <- paste("taskset -c 0-11 Rscript", file.path(script_dir, mtd), siv)
        
        log_file <- file.path(log_dir, paste0(mtd, ".i.SI", round(siv, 2), ".log"))
        cmdq <- paste(cmd, ">", log_file, "2>&1")

        exit_code <- system(cmdq)
        if(exit_code != 0) {
            err_msg <- paste0("!!! [Worker] Error in method ", mtd, " at SI ", round(siv, 2), "\n")
            cat(err_msg)
        }
    }

    return(paste0("SI ", siv, " Completed."))
}

stopCluster(cl)
gc()
cat("\n Islet Done. \n")
cat("\n=== All simulations completed! ===\n") 
cat("\n Elapsed time: ", difftime(Sys.time(), t_start, units = "hours"), "\n")

