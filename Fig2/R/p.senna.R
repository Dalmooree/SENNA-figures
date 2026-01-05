suppressPackageStartupMessages(
    suppressWarnings( {
        library(SENNA)
        library(Seurat)
        library(dplyr)
    }))

args <- commandArgs(trailingOnly = TRUE)
si <- as.numeric(args[1])
sit <- round(si, 2)
dpath <- "data/sim.generated"
rpath <- "data/bhmk.res/p/SENNA"
spath <- "data/Thymus9/outs"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

prog_count <- read.csv(file.path(dpath, "sim.prog", paste0("sim_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
prog_perm <- read.csv(file.path(dpath, "sim.prog", paste0("perm_prog_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.prog", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
gs <- g$gs; gc <- g$gc
pval <- 0.05






# SENNA regio -----------------------------------------------

knot <- read.csv("data/knots/svgp.csv", header = TRUE)

seurat <- Load10X_Spatial(
  data.dir = spath,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Visium",
  filter.matrix = TRUE,
  image = NULL)

if(min(seurat$nCount_Spatial) == 0) seurat <- subset(seurat, nCount_Spatial > 0)

t0 <- Sys.time()

surt_prog <- CreateSeuratObject(
    counts = as.sparse(prog_count),
    assay = "Spatial",
    meta.data = seurat@meta.data)
surt_prog@images <- seurat@images

surt_prog@meta.data$nCount_Spatial <- colSums(prog_count)
surt_prog@meta.data$nFeature_Spatial <- colSums(prog_count > 0)
surt_prog <- NormalizeData(surt_prog, verbose = FALSE)
surt_prog <- ScaleData(surt_prog, 
                        features = c(gs, gc),
                        verbose = FALSE)


sensim <- SENNA_Visium(surt_prog,
                       assay = "Spatial",
                       all_genes = TRUE,
                       slice_name = "Visium")

rm(surt_prog); gc()
                       
sensim <- FullCurve(senna = sensim,
                    knot_df = knot,
                    type = "spline")
sensim <- GetCurveParam(sensim)

t1 <- Sys.time()

sensim <- ProgSVGs(sensim,
                   FDR_level = pval,
                   grad_cutoff = 0.1)
sensim <- c(sensim@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]],
            sensim@Gene[["P.SVGs"]][["Variable_gene"]][["negative"]])

t2 <- Sys.time()
  
surt_perm <- CreateSeuratObject(
  counts = as.sparse(prog_perm),
  assay = "Spatial",
  meta.data = seurat@meta.data)
surt_perm@images <- seurat@images

rm(seurat); rm(surt_prog); gc()

surt_perm@meta.data$nCount_Spatial <- colSums(prog_perm)
surt_perm@meta.data$nFeature_Spatial <- colSums(prog_perm > 0)
DefaultAssay(surt_perm) <- "Spatial"
surt_perm <- NormalizeData(surt_perm, verbose = FALSE)
surt_perm <- ScaleData(surt_perm, 
                       features = c(gs, gc),
                       verbose = FALSE)
senprm <- SENNA_Visium(surt_perm,
                       assay = "Spatial",
                       all_genes = TRUE,
                       slice_name = "Visium")

rm(surt_perm); gc()

senprm <- FullCurve(senna = senprm,
                    knot_df = knot,
                    type = "spline")
senprm <- GetCurveParam(senprm)
senprm <- ProgSVGs(senprm,
                   FDR_level = pval,
                   grad_cutoff = 0.1)
senprm <- c(senprm@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]],
            senprm@Gene[["P.SVGs"]][["Variable_gene"]][["negative"]])

# Save results ------------------------------------------------------------

res <- tibble(
    ID = c("Power_sim", "FPR_sim", "Power_prm", "FPR_prm"),
    SENNA = c(length(intersect(gs, sensim)) / length(gs),
                length(intersect(gc, sensim)) / length(gc),
                length(intersect(gs, senprm)) / length(gs),
                length(intersect(gc, senprm)) / length(gc)),
    SI = rep(sit, 4),
    total_time = rep(as.numeric(difftime(t2, t0, units = "secs")), 4),
    mapping_time = rep(as.numeric(difftime(t1, t0, units = "secs")), 4)
)

saveRDS(res, file = paste0(rpath, "/res_senna_SI", sit, ".rds"))