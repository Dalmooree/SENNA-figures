suppressPackageStartupMessages(
    suppressWarnings( {
        library(SPARK)
        library(dplyr)
    }))

args <- commandArgs(trailingOnly = TRUE)
si <- as.numeric(args[1])
sit <- round(si, 2)
dpath <- "data/sim.generated"
rpath <- "data/bhmk.res/p/SPARKX"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

prog_count <- read.csv(file.path(dpath, "sim.prog", paste0("sim_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
prog_perm <- read.csv(file.path(dpath, "sim.prog", paste0("perm_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.prog", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
gs <- g$gs; gc <- g$gc

sparkarg <- read.csv("data/sim.generated/coord.csv", header = TRUE, row.names = 1)
pval <- 0.05



# SPARKX regio -----------------------------------------------

sparkarg <- as.matrix(sparkarg)
prog_count <- as.matrix(prog_count)
prog_perm <- as.matrix(prog_perm)

t0 <- Sys.time()

spksim <- sparkx(prog_count,
                 sparkarg,
                 option ="mixture", 
                 verbose = FALSE)
spksim <- rownames(dplyr::filter(spksim$res_mtest, 
                                 adjustedPval <= pval))

t1 <- Sys.time()

spkprm <- sparkx(prog_perm,
                 sparkarg,
                 option ="mixture", 
                 verbose = FALSE)
spkprm <- rownames(dplyr::filter(spkprm$res_mtest, 
                                 adjustedPval <= pval))


# Save results ------------------------------------------------------------
res <- tibble(
    ID = c("Power_sim", "FPR_sim", "Power_prm", "FPR_prm"),
    SPARK_X = c(length(intersect(gs, spksim)) / length(gs),
                length(intersect(gc, spksim)) / length(gc),
                length(intersect(gs, spkprm)) / length(gs),
                length(intersect(gc, spkprm)) / length(gc)),
    SI = rep(sit, 4),
    total_time = rep(as.numeric(difftime(t1, t0, units = "secs")), 4),
    mapping_time = rep(NA, 4)
)

saveRDS(res, 
        file = paste0(rpath, "/res_sparkx_SI", sit, ".rds"))