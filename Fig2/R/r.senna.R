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
rpath <- "data/bhmk.res/r/SENNA"
spath <- "data/Thymus9/outs"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

regio_count <- read.csv(file.path(dpath, "sim.regio", paste0("sim_regio_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
regio_perm <- read.csv(file.path(dpath, "sim.regio", paste0("perm_regio_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.regio", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
gs <- g$gs; gc <- g$gc
pval <- 0.05






# SENNA regio -----------------------------------------------

knot <- read.csv("data/knots/svgr.csv", header = TRUE)

seurat <- Load10X_Spatial(
  data.dir = spath,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Visium",
  filter.matrix = TRUE,
  image = NULL)

if(min(seurat$nCount_Spatial) == 0) seurat <- subset(seurat, nCount_Spatial > 0)

t0 <- Sys.time()

surt_regio <- CreateSeuratObject(
    counts = as.sparse(regio_count),
    assay = "Spatial",
    meta.data = seurat@meta.data)
surt_regio@images <- seurat@images

surt_regio@meta.data$nCount_Spatial <- colSums(regio_count)
surt_regio@meta.data$nFeature_Spatial <- colSums(regio_count > 0)
surt_regio <- NormalizeData(surt_regio, verbose = FALSE)
surt_regio <- ScaleData(surt_regio, 
                        features = c(gs, gc),
                        verbose = FALSE)


sensim <- SENNA_Visium(surt_regio,
                       assay = "Spatial",
                       all_genes = TRUE,
                       slice_name = "Visium")

rm(surt_regio); gc()
                       
sensim <- FullCurve(senna = sensim,
                    knot_df = knot,
                    type = "spline")
sensim <- GetCurveParam(sensim)
sensim <- TissueRegionation(sensim)

t1 <- Sys.time()

sensim <- RegionSVGs(sensim,
                     FDR_level = pval,
                     grad_cutoff = 0.1)
sensim <- c(sensim@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
            sensim@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])

t2 <- Sys.time()
  
surt_perm <- CreateSeuratObject(
  counts = as.sparse(regio_perm),
  assay = "Spatial",
  meta.data = seurat@meta.data)
surt_perm@images <- seurat@images

rm(seurat); rm(surt_regio); gc()

surt_perm@meta.data$nCount_Spatial <- colSums(regio_perm)
surt_perm@meta.data$nFeature_Spatial <- colSums(regio_perm > 0)
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
senprm <- TissueRegionation(senprm)
senprm <- RegionSVGs(senprm,
                     FDR_level = pval,
                     grad_cutoff = 0.1)
senprm <- c(senprm@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
            senprm@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])


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