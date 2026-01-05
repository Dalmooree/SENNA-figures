suppressPackageStartupMessages(
    suppressWarnings( {
        library(Seurat)
        library(dplyr)
    }))

args <- commandArgs(trailingOnly = TRUE)
si <- as.numeric(args[1])
sit <- round(si, 2)
dpath <- "data/sim.generated"
rpath <- "data/bhmk.res/i/SEURAT"
spath <- "data/Thymus9/outs"
if(!dir.exists(rpath)) dir.create(rpath, recursive = TRUE)

isl_count <- read.csv(file.path(dpath, "sim.isl", paste0("sim_isl_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
isl_perm <- read.csv(file.path(dpath, "sim.isl", paste0("perm_isl_SI", sit, ".csv")), row.names = 1, header = TRUE, check.names = FALSE)
g <- read.csv(file.path(dpath, "sim.isl", paste0("GT_SI", sit, ".csv")), header = TRUE, check.names = FALSE)
gs <- g$gs; gc <- g$gc
pval <- 0.05


# Seurat Moran's I -----------------------------------------------

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

morsim <- Seurat::RunMoransI(
  data = surt_isl[["Spatial"]]$scale.data,
  pos = GetTissueCoordinates(object = surt_isl[["Visium"]])[,1:2],
  verbose = FALSE)
morsim$p.value <- p.adjust(morsim$p.value, method = "BH")
morsim <- rownames(dplyr::filter(morsim, p.value <= pval))

t1 <- Sys.time()


surt_perm <- CreateSeuratObject(
    counts = as.sparse(isl_perm),
    assay = "Spatial",
    meta.data = seurat@meta.data)
surt_perm@images <- seurat@images

surt_perm@meta.data$nCount_Spatial <- colSums(isl_perm)
surt_perm@meta.data$nFeature_Spatial <- colSums(isl_perm > 0)
DefaultAssay(surt_perm) <- "Spatial"
surt_perm <- NormalizeData(surt_perm, verbose = FALSE)
surt_perm <- ScaleData(surt_perm, 
                       features = c(gs, gc),
                       verbose = FALSE)

rm(seurat); gc()
  
rm(surt_isl); gc()
  
morprm <- Seurat::RunMoransI(
  data = surt_perm[["Spatial"]]$scale.data,
  pos = GetTissueCoordinates(object = surt_perm[["Visium"]])[,1:2],
  verbose = FALSE)
morprm$p.value <- p.adjust(morprm$p.value, method = "BH")
morprm <- rownames(dplyr::filter(morprm, p.value <= pval))

rm(surt_perm); gc()


# Save results ------------------------------------------------------------

res <- tibble(
    ID = c("Power_sim", "FPR_sim", "Power_prm", "FPR_prm"),
    SEURAT = c(length(intersect(gs, morsim)) / length(gs),
                length(intersect(gc, morsim)) / length(gc),
                length(intersect(gs, morprm)) / length(gs),
                length(intersect(gc, morprm)) / length(gc)),
    SI = rep(sit, 4),
    total_time = rep(as.numeric(difftime(t1, t0, units = "secs")), 4),
    mapping_time = rep(NA, 4))

saveRDS(res, 
        file = paste0(rpath, "/res_seurat_SI", sit, ".rds"))