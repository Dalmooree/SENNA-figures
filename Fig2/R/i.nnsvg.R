suppressPackageStartupMessages(
    suppressWarnings( {
        library(SpatialExperiment)
        library(scran)
        library(nnSVG)
        library(dplyr)
    }))

args <- commandArgs(trailingOnly = TRUE)
si <- as.numeric(args[1])
sit <- round(si, 2)
dpath <- "data/sim.generated"
rpath <- "data/bhmk.res/i/NNSVG"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

isl_count <- read.csv(file.path(dpath, "sim.isl", paste0("sim_isl_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
isl_perm <- read.csv(file.path(dpath, "sim.isl", paste0("perm_isl_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.isl", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
gs <- g$gs; gc <- g$gc
pval <- 0.05

# NNSVG isl -----------------------------------------------
isl_count <- as.matrix(isl_count)
isl_perm <- as.matrix(isl_perm)
coords <- read.csv("data/sim.generated/coord.csv", header = TRUE, row.names = 1)
coords <- as.matrix(coords)
colnames(coords) <- c("x", "y")

t0 <- Sys.time()

spe <- SpatialExperiment(
    assays = list(counts = isl_count),
    spatialCoords = coords
)
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
spe <- nnSVG(spe,
             assay_name = "logcounts",
             n_neighbors = 10,
             n_threads = 1,
             verbose = FALSE)


svgsim <- rownames(rowData(spe))[which(rowData(spe)$padj <= pval)]

t1 <- Sys.time()

spe <- SpatialExperiment(
    assays = list(counts = isl_perm),
    spatialCoords = coords
)
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
spe <- nnSVG(spe,
             assay_name = "logcounts",
             n_neighbors = 10,
             n_threads = 1,
             verbose = FALSE)

svgprm <- rownames(rowData(spe))[which(rowData(spe)$padj <= pval)]

# Save results ------------------------------------------------------------
res <- tibble(
    ID = c("Power_sim", "FPR_sim", "Power_prm", "FPR_prm"),
    NNSVG = c(length(intersect(gs, svgsim)) / length(gs),
              length(intersect(gc, svgsim)) / length(gc),
              length(intersect(gs, svgprm)) / length(gs),
              length(intersect(gc, svgprm)) / length(gc)),
    SI = rep(sit, 4),
    total_time = rep(as.numeric(difftime(t1, t0, units = "secs")), 4),
    mapping_time = rep(NA, 4)
)

saveRDS(res, 
        file = paste0(rpath, "/res_nnsvg_SI", sit, ".rds"))
