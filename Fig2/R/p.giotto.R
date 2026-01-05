suppressPackageStartupMessages(
    suppressWarnings( {
        library(Giotto)
        library(Seurat)
        library(dplyr)
    }))

args <- commandArgs(trailingOnly = TRUE)
si <- as.numeric(args[1])
sit <- round(si, 2)
dpath <- "data/sim.generated"
rpath <- "data/bhmk.res/p/GIOTTO"
spath <- "data/Thymus9/outs"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

prog_count <- read.csv(file.path(dpath, "sim.prog", paste0("sim_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
prog_perm <- read.csv(file.path(dpath, "sim.prog", paste0("perm_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.prog", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
coords <- read.csv("data/sim.generated/coord.csv", header = TRUE, row.names = 1)
gs <- g$gs; gc <- g$gc
pval <- 0.05


# Giotto kmeans -----------------------------------------------

t0 <- Sys.time()
suppressWarnings({
  giksim <- createGiottoObject(expression = as.matrix(prog_count),
    spatial_locs = coords,
    instructions = createGiottoInstructions(show_plot = FALSE, save_plot = FALSE))
  giksim <- normalizeGiotto(gobject = giksim)
  giksim <- createSpatialNetwork(gobject = giksim, minimum_k = 0)
  giksim <- binSpect(giksim, bin_method = 'kmeans', verbose = FALSE)
  giksim <- dplyr::filter(giksim, adj.p.value <= pval)[["feats"]]
})

t1 <- Sys.time()

suppressWarnings({
    gikprm <- createGiottoObject(expression = as.matrix(prog_perm),
      spatial_locs = coords,
      instructions = createGiottoInstructions(show_plot = FALSE, save_plot = FALSE))
    gikprm <- normalizeGiotto(gobject = gikprm)
    gikprm <- createSpatialNetwork(gobject = gikprm, minimum_k = 0)
    gikprm <- binSpect(gikprm, bin_method = 'kmeans', verbose = FALSE)
    gikprm <- dplyr::filter(gikprm  , adj.p.value <= pval)[["feats"]]
})

# Save results -----------------------------------------------

res <- tibble(
    ID = c("Power_sim", "FPR_sim", "Power_prm", "FPR_prm"),
    GIOTTO_kmeans = c(length(intersect(gs, giksim)) / length(gs),
                      length(intersect(gc, giksim)) / length(gc),
                      length(intersect(gs, gikprm)) / length(gs),
                      length(intersect(gc, gikprm)) / length(gc)),
    SI = rep(sit, 4),
    total_time = rep(as.numeric(difftime(t1, t0, units = "secs")), 4),
    mapping_time = rep(NA, 4))

saveRDS(res, 
        file = paste0(rpath, "/res_giottok_SI", sit, ".rds"))



# Giotto rank -----------------------------------------------

t0 <- Sys.time()
suppressWarnings({
  girsim <- createGiottoObject(expression = as.matrix(prog_count),
    spatial_locs = coords,
    instructions = createGiottoInstructions(show_plot = FALSE, save_plot = FALSE))
  girsim <- normalizeGiotto(gobject = girsim)
  girsim <- createSpatialNetwork(gobject = girsim, minimum_k = 0)
  girsim <- binSpect(girsim, bin_method = 'rank', verbose = FALSE)
  girsim <- dplyr::filter(girsim, adj.p.value <= pval)[["feats"]]
})

t1 <- Sys.time()

suppressWarnings({
    girprm <- createGiottoObject(expression = as.matrix(prog_perm),
      spatial_locs = coords,
      instructions = createGiottoInstructions(show_plot = FALSE, save_plot = FALSE))
    girprm <- normalizeGiotto(gobject = girprm)
    girprm <- createSpatialNetwork(gobject = girprm, minimum_k = 0)
    girprm <- binSpect(girprm, bin_method = 'rank', verbose = FALSE)
    girprm <- dplyr::filter(girprm  , adj.p.value <= pval)[["feats"]]
})

# Save results -----------------------------------------------

res <- tibble(
    ID = c("Power_sim", "FPR_sim", "Power_prm", "FPR_prm"),
    GIOTTO_rank = c(length(intersect(gs, girsim)) / length(gs),
                      length(intersect(gc, girsim)) / length(gc),
                      length(intersect(gs, girprm)) / length(gs),
                      length(intersect(gc, girprm)) / length(gc)),
    SI = rep(sit, 4),
    total_time = rep(as.numeric(difftime(t1, t0, units = "secs")), 4),
    mapping_time = rep(NA, 4))

saveRDS(res, 
        file = paste0(rpath, "/res_giottor_SI", sit, ".rds"))
