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
rpath <- "data/bhmk.res/i/SENNA"
spath <- "data/Thymus9/outs"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

isl_count <- read.csv(file.path(dpath, "sim.isl", paste0("sim_isl_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
isl_perm <- read.csv(file.path(dpath, "sim.isl", paste0("perm_isl_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.isl", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
gs <- g$gs; gc <- g$gc
pval <- 0.05






# SENNA regio -----------------------------------------------

knot <- read.csv("data/knots/svgi.csv", header = TRUE)

seurat <- Load10X_Spatial(
  data.dir = spath,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Visium",
  filter.matrix = TRUE,
  image = NULL)

if(min(seurat$nCount_Spatial) == 0) seurat <- subset(seurat, nCount_Spatial > 0)

t0 <- Sys.time()

surt_isl <- CreateSeuratObject(
    counts = as.sparse(isl_count),
    assay = "Spatial",
    meta.data = seurat@meta.data)
surt_isl@images <- seurat@images

surt_isl@meta.data$nCount_Spatial <- colSums(isl_count)
surt_isl@meta.data$nFeature_Spatial <- colSums(isl_count > 0)
surt_isl <- NormalizeData(surt_isl, verbose = FALSE)
surt_isl <- ScaleData(surt_isl, 
                        features = c(gs, gc),
                        verbose = FALSE)


sensim <- SENNA_Visium(surt_isl,
                       assay = "Spatial",
                       all_genes = TRUE,
                       slice_name = "Visium")

rm(surt_isl); gc()
                       
sensim <- TrimmedCurve(senna = sensim,
                       knot_df = knot,
                       type = "islet")
sensim <- GetCurveParam(sensim)
sensim <- TissueRegionation(sensim)

t1 <- Sys.time()

sensim <- RegionSVGs(sensim,
                     direction = -1,
                     FDR_level = pval,
                     grad_cutoff = 0.1)
sensim <- c(sensim@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
            sensim@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])   
t2 <- Sys.time()
  
surt_perm <- CreateSeuratObject(
  counts = as.sparse(isl_perm),
  assay = "Spatial",
  meta.data = seurat@meta.data)
surt_perm@images <- seurat@images

rm(seurat); rm(surt_isl); gc()

surt_perm@meta.data$nCount_Spatial <- colSums(isl_perm)
surt_perm@meta.data$nFeature_Spatial <- colSums(isl_perm > 0)
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

senprm <- TrimmedCurve(senna = senprm,
                       knot_df = knot,
                       type = "islet")
senprm <- GetCurveParam(senprm)
senprm <- TissueRegionation(senprm)
senprm <- RegionSVGs(senprm,
                     direction = -1,
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