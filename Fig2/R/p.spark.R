suppressPackageStartupMessages(
    suppressWarnings( {
        library(SPARK)
        library(dplyr)
    }))

args <- commandArgs(trailingOnly = TRUE)
si <- as.numeric(args[1])
sit <- round(si, 2)
dpath <- "data/sim.generated"
rpath <- "data/bhmk.res/p/SPARK"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

prog_count <- read.csv(file.path(dpath, "sim.prog", paste0("sim_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
prog_perm <- read.csv(file.path(dpath, "sim.prog", paste0("perm_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.prog", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
gs <- g$gs; gc <- g$gc

sparkarg <- read.csv("data/sim.generated/coord.csv", header = TRUE, row.names = 1)

pval <- 0.05


# SPARK regio -----------------------------------------------

t0 <- Sys.time()
spsim <- CreateSPARKObject(counts = prog_count,
                            location = sparkarg,
                            percentage = 0.1)
spsim@lib_size <- apply(spsim@counts, 2, sum)
spsim <- spark.vc(spsim, 
                  covariates = NULL, 
                  lib_size = spsim@lib_size, 
                  num_core = 1L,
                  verbose = F)
spsim <- spark.test(spsim, 
                    check_positive = T, 
                    verbose = F)
spsim <- spsim@res_mtest
spsim <- rownames(dplyr::filter(spsim, adjusted_pvalue <= pval))

t1 <- Sys.time()            


spprm <- CreateSPARKObject(counts = prog_perm,
                            location = sparkarg,
                            percentage = 0.1)
spprm@lib_size <- apply(spprm@counts, 2, sum)
spprm <- spark.vc(spprm, 
                  covariates = NULL, 
                  lib_size = spprm@lib_size, 
                  num_core = 1L,
                  verbose = F)
spprm <- spark.test(spprm, 
                    check_positive = T, 
                    verbose = F)
spprm <- spprm@res_mtest
spprm <- rownames(dplyr::filter(spprm, adjusted_pvalue <= pval))



# Save results ------------------------------------------------------------

res <- tibble(
    ID = c("Power_sim", "FPR_sim", "Power_prm", "FPR_prm"),
    SPARK = c(length(intersect(gs, spsim)) / length(gs),
              length(intersect(gc, spsim)) / length(gc),
              length(intersect(gs, spprm)) / length(gs),
              length(intersect(gc, spprm)) / length(gc)),
    SI = rep(sit, 4),
    total_time = rep(as.numeric(difftime(t1, t0, units = "secs")), 4),
    mapping_time = rep(NA, 4)
)

saveRDS(res, 
        file = paste0(rpath, "/res_spark_SI", sit, ".rds"))

