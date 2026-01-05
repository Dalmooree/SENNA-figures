suppressPackageStartupMessages(
  suppressWarnings({
    library(Seurat)
    library(SENNA)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(tidyr)
  })
)

# path
path <- getwd()


# Spatial ATAC-RNA co-profiling ----
rescalecp <- function(senna){
  coord_mat <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
  fx <- senna@CurveAxis[["fun"]][["x.coef"]]
  fy <- senna@CurveAxis[["fun"]][["y.coef"]]
  knots <- senna@CurveAxis[["knots"]]
  kn <- senna@CurveAxis[["fun"]][["t"]]
  
  if("trimmed" %in% senna@CurveAxis[["type"]]) {
    ll <- lapply(kn[-length(kn)],
                 function(k){
                   qf <- function(x){
                     sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                            (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                   }
                   return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                 })
    ll <- c(0, unlist(ll))
    cll <- cumsum(ll)
    tv <- coord_mat[["t"]]
    tprime <- 
      ifelse(is.nan(tv),
             tv,
             {
               fidx <- floor(tv)
               fidx[fidx == 0] <- 1L
               fidx[fidx > length(cll)] <- length(cll)
               
               (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                 cll[length(cll)]
             })
    senna@Coord[["Spatial"]][["tprime"]] <- NA_real_
    senna@Coord[["Spatial"]][["tprime"]][!is.na(senna@Coord[["Spatial"]]$t)] <- tprime
  } else {
    kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
    tv <- coord_mat[["t"]]
    tmin <- min(tv); tmax <- max(tv)
    llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                    ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
    lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                    ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
    ll <- lapply(2:(mk-1),
                 function(k){
                   qf <- function(x){
                     sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                            (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                   }
                   return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                 })
    ll <- c(0, llw, unlist(ll), lhi)
    cll <- cumsum(ll)
    tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
    tv <- ifelse(tv > (max(kn) - 1L), 
                 (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                 tv)
    tv <- tv + 1L
    tprime <- {
      fidx <- floor(tv)
      fidx[fidx == 0] <- 1L
      fidx[fidx >= length(cll)] <- length(cll) - 1L
      (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
        cll[length(cll)]
    }
    senna@Coord[["Spatial"]][["tprime"]] <- tprime
  }
  
  return(senna)
}

##

# Load Data
dpath <- file.path(getwd(), "dataset", "UCSCcellbrowser_P22mousebrain_spatial_RNA_ATAC")
z <- readRDS(file.path(dpath, "P22mousebrain_spatial_RNA_ATAC.rds"))

# SENNA ----
feats <- intersect(Features(z, assay = "SCT"), Features(z, assay = "ATAC"))
DefaultAssay(z) <- "Spatial"
z <- NormalizeData(z, normalization.method = "LogNormalize", verbose = FALSE)
z <- ScaleData(z, verbose = FALSE)
VariableFeatures(z) <- feats
Idents(z) <- z$RNA_clusters
sen_rna <- SENNA_Visium(z, slice_name = "slice1", annotation = TRUE)

z@assays$ATAC@scale.data <- as.matrix(z@assays$ATAC@data)
Idents(z) <- z$ATAC_clusters
VariableFeatures(z, assay = "ATAC") <- feats
sen_atac <- SENNA_Visium(z, assay = "ATAC", slice_name = "slice1", annotation = TRUE)

tmp <- sen_rna@Coord$Spatial
sen_rna@Coord$Spatial$X1 <- 
  sen_atac@Coord$Spatial$X1 <- tmp$X2
sen_rna@Coord$Spatial$X2 <- 
  sen_atac@Coord$Spatial$X2 <- 1 - tmp$X1

TissueFOV(sen_rna)
TissueFOV(sen_atac)

rm("z", "tmp"); gc()

### pre-processing

#AppDat(sen_rna, reference_value = "Annotation", colorset = Polychrome::alphabet.colors(11))
#knot_picker()

knots <- read.csv(file.path(dpath, "cortex.csv"))
sen_rna <- TrimmedCurve(sen_rna, knots, type = "spline")
sen_atac <- TrimmedCurve(sen_atac, knots, type = "spline")

sen_rna <- GetCurveParam(sen_rna)
sen_atac <- GetCurveParam(sen_atac)

griddat <- crvtrjry(sen_rna)
ScanInterval(sen_rna, interval = .13, annotation = TRUE, dot_size = 0.6) + 
  scale_color_manual(values = unname(Polychrome::alphabet.colors(11))) +
  geom_point(aes(X1, X2), data = griddat, color = "black", size = .1) +
  geom_point(aes(X1, X2), data = knots, color = "black") +
  ggrepel::geom_text_repel(aes(X1, X2), label = 1:nrow(knots), data = knots) + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())


ppath <- file.path(dpath, "Fig")
ggsave(file.path(ppath, "intervals.png"), width = 4, height = 4, units = "in")


sen_rna <- ProgSVGs(sen_rna,
                    FDR_level = 0.01,
                    grad_cutoff = .1,
                    interval = .13)
sen_atac <- ProgSVGs(sen_atac,
                     FDR_level = 0.01,
                     grad_cutoff = .01,
                     interval = .13)
#saveRDS(sen_rna, file.path(dpath, "sen_rna_fitted.rds"))
#saveRDS(sen_atac, file.path(dpath, "sen_atac_fitted.rds"))

# Output----
dpath <- file.path(getwd(), "dataset", "UCSCcellbrowser_P22mousebrain_spatial_RNA_ATAC")
ppath <- file.path(dpath, "Fig")

sen_rna <- readRDS(file.path(dpath, "sen_rna_fitted.rds"))
sen_atac <- readRDS(file.path(dpath, "sen_atac_fitted.rds"))


## Fig.5a
sen_rna@Gene[["Reference"]][["Annotation"]] <- factor(
  sen_rna@Gene[["Reference"]][["Annotation"]], 
  c("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
)
sen_atac@Gene[["Reference"]][["Annotation"]] <- factor(
  sen_atac@Gene[["Reference"]][["Annotation"]], 
  c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13")
)

p5a1 <- TissueFOV(sen_rna,
                  dot_size = .3,
                  dot_alpha = .8,
                  colors = unname(Polychrome::sky.colors(
                    length(unique(sen_rna@Gene[["Reference"]][["Annotation"]]))
                  ))) +
  labs(title = "Cluster_RNA",
       color = "Cluster") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        title = element_text(size = 14),
        legend.margin = margin(l = -3, r = 3),
        legend.box.margin = margin(l = -3, r = 3),
        legend.spacing.y = unit(0, units = "pt"),
        legend.box.just = "left",
        legend.key.height = unit(14, "pt")) +
  guides(color = guide_legend(override.aes = list(size = 2,
                                                  alpha = 1),
                              byrow = TRUE)); p5a1


p5a2 <- TissueFOV(sen_atac,
                  dot_size = .3,
                  dot_alpha = .8,
                  colors = unname(Polychrome::kelly.colors(
                    length(unique(sen_atac@Gene[["Reference"]][["Annotation"]]))
                  ))) +
  labs(title = "Cluster_ATAC_GAS",
       color = "Cluster") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        title = element_text(size = 14),
        legend.margin = margin(l = -3),
        legend.box.margin = margin(l = -3),
        legend.spacing.y = unit(-2, units = "pt"),
        legend.box.just = "left",
        legend.key.height = unit(14, "pt")) +
  guides(color = guide_legend(override.aes = list(size = 2,
                                                  alpha = 1),
                              byrow = TRUE)); p5a2

p5a <- p5a2 + p5a
ggsave(plot = p5a,
       filename = file.path(path, "SENNA", "Fig", "5_multiome", "ATAC", "cls.tif"),
       width = 9, height = 4 , dpi = 600)




feats <- c("Mef2c", "Neurod6")

datr <- sen_rna@Coord[["Spatial"]] %>%
  mutate(aph = is.finite(.$distance) & abs(.$distance) < 0.13)
datr <- cbind(datr, sen_rna@Gene[["Spatial"]][, feats])
data <- sen_atac@Coord[["Spatial"]] %>%
  mutate(aph = is.finite(.$distance) & abs(.$distance) < 0.13)
data <- cbind(data, sen_atac@Gene[["Spatial"]][, feats])

resr <- filter(sen_rna@Gene[["P.SVGs"]][["Report"]], Gene %in% feats)
resa <- filter(sen_atac@Gene[["P.SVGs"]][["Report"]], Gene %in% feats)



pmr <- ggplot() +
  geom_point(aes(X1, X2, color = .data[[feats[1]]], alpha = aph), data = datr,
             size = .6) +
  viridis::scale_color_viridis() + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6)) +
  labs(title = expression(italic("Mef2c")),
       color = "RNA") +
  geom_path(aes(X1, X2),
            data = crvtrjry(sen_rna),
            color = "#000000",
            size = 1) +
  theme_light() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        title = element_text(size = 14),
        legend.margin = margin(l = -3),
        legend.box.margin = margin(l = -3),
        legend.box.just = "left",
        aspect.ratio = 1) +
  guides(alpha = "none",
         color = guide_colorbar(barwidth = 1,
                                barheight = 6)) ; pmr

pnr <- ggplot() +
  geom_point(aes(X1, X2, color = .data[[feats[2]]], alpha = aph), data = datr,
             size = .6) +
  viridis::scale_color_viridis() + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6)) +
  labs(title = expression(italic("Neurod6")),
       color = "RNA") +
  geom_path(aes(X1, X2),
            data = crvtrjry(sen_rna),
            color = "#000000",
            size = 1) +
  theme_light() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        title = element_text(size = 14),
        legend.margin = margin(l = -3),
        legend.box.margin = margin(l = -3),
        legend.box.just = "left",
        aspect.ratio = 1) + 
  guides(alpha = "none",
         color = guide_colorbar(barwidth = 1,
                                barheight = 6)) ; pnr


pb <- (pmr / pnr)
ggsave(filename = file.path(path, "SENNA", "Fig", "5_multiome", "ATAC", "rna.tif"), 
       plot = pb, width = 4.58, height = 8.52, dpi = 600)



pma <- ggplot() +
  geom_point(aes(X1, X2, color = .data[[feats[1]]], alpha = aph), data = data,
             size = .6) +
  viridis::scale_color_viridis(option = "turbo") + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6)) +
  labs(title = expression(italic("Mef2c")),
       color = "GAS") +
  geom_path(aes(X1, X2),
            data = crvtrjry(sen_atac),
            color = "#000000",
            size = 1) +
  theme_light() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        title = element_text(size = 14),
        legend.margin = margin(l = -3),
        legend.box.margin = margin(l = -3),
        legend.box.just = "left",
        aspect.ratio = 1) +
  guides(alpha = "none",
         color = guide_colorbar(barwidth = 1,
                                barheight = 6)) ; pma

pna <- ggplot() +
  geom_point(aes(X1, X2, color = .data[[feats[2]]], alpha = aph), data = data,
             size = .6) +
  viridis::scale_color_viridis(option = "turbo") + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6)) +
  labs(title = expression(italic("Neurod6")),
       color = "GAS") +
  geom_path(aes(X1, X2),
            data = crvtrjry(sen_atac),
            color = "#000000",
            size = 1) +
  theme_light() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        title = element_text(size = 14),
        legend.margin = margin(l = -3),
        legend.box.margin = margin(l = -3),
        legend.box.just = "left",
        aspect.ratio = 1) + 
  guides(alpha = "none",
         color = guide_colorbar(barwidth = 1,
                                barheight = 6)) ; pna


pc <- (pma / pna)
ggsave(filename = file.path(path, "SENNA", "Fig", "5_multiome", "ATAC", "gas.tif"), 
       plot = pc, width = 4.58, height = 8.52, dpi = 600)











# CODEX preprocessing
cod <- LoadAkoya(
  filename = paste0(path,
                    "dataset/publicCODEX/HBM754_WKLP_262/LN7910_20_008_11022020_reg001_compensated.csv"),
  type = "processor", 
  fov = "HBM754.WKLP.262")

cod <- NormalizeData(object = cod, 
                     normalization.method = "CLR", 
                     margin = 2,
                     verbose = FALSE)
cod <- ScaleData(cod, verbose = FALSE)
VariableFeatures(cod) <- rownames(cod)
cod <- RunPCA(object = cod, npcs = 10, verbose = FALSE)
cod <- RunUMAP(object = cod, dims = 1:10, verbose = FALSE)
cod <- FindNeighbors(object = cod, dims = 1:10, verbose = FALSE)
cod <- FindClusters(object = cod, verbose = FALSE, resolution = 0.4, n.start = 1)

ImageDimPlot(cod, 
             cols = DiscretePalette(
               length(unique(cod$seurat_clusters)), 
               palette = "parade", 
               shuffle = FALSE)) +
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(fill = "Cluster")

ggsave(
  plot = p5a1,
  filename = paste0(
    path,
    "SENNA/Fig/5_codex/seuratdimplot1.tif"
    ), 
  height = 8 * .6, width = 6 * .6, dpi = 600)



#saveRDS(cod, paste0(path, "dataset/rds/codex_seurat.rds"))


# SENNA
cod <- readRDS(paste0(path, "dataset/rds/codex_seurat.rds"))

sen <- SENNA_CODEX(cod,
                   assay = "Akoya",
                   fov_name = "HBM754.WKLP.262",
                   annotation = TRUE)

#AppDat(sen, reference_value = "Annotation", colorset = colorset)
#knot_picker()

ck <- read.csv(
  paste0(path, "dataset/knots/codex_knots.csv"), 
  header = TRUE)

sen <- FullCurve(sen, ck, "spline")



sen <- readRDS(paste0(path, "dataset/rds/codex_senna.rds"))

plt_factor <- (max(sen@Coord$Spatial$X2) - min(sen@Coord$Spatial$X2)) / 
  (max(sen@Coord$Spatial$X1) - min(sen@Coord$Spatial$X1))

colorset <- DiscretePalette(
  length(unique(sen@Gene$Reference$Annotation)), 
  palette = "parade", 
  shuffle = FALSE)

p5d <- ShowCurve(sen,
                 color_reference = "Annotation",
                 colors = colorset,
                 line_color = "#000000",
                 knots_color = "#000000",
                 bg_dot_size = .1,
                 bg_dot_alpha = .3,
                 knots_size = 0,
                 line_size = 1,
                 order_label = FALSE) +
  theme_light() +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +
  labs(color = "Cluster") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -5),
        legend.box.just = "left",
        legend.key.height = unit(15, "pt"),
        aspect.ratio = plt_factor); p5d


ggsave(plot = p5d,
       filename = paste0(
         path,
         "SENNA/Fig/5_multiome/LN/lncls.tif"
       ), 
       width = 4.2, height = 6.2, dpi = 600)



sen <- GetCurveParam(sen)
sen <- TissueRegionation(sen)

p5e <- ShowCSDistance(sen,
                      dot_size = .4,
                      dot_alpha = .3) +
  scale_color_gradient2(high = "#32a02c",
                        mid = "#dedede",
                        low = "#6a3d9a",
                        breaks = c(min(sen@Coord$Spatial$distance) * .8, 
                                   0, 
                                   max(sen@Coord$Spatial$distance) * .8),
                        labels = c("-", 0, "+")) +
  labs(color = "CSD") +
  theme_light() +
  guides(color = guide_colorbar(barwidth = 1,
                                barheight = 6)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -5),
        legend.box.just = "left",
        aspect.ratio = plt_factor); p5e

ggsave(
  plot = p5e,
  filename = paste0(
    path,
    "SENNA/Fig/5_multiome/LN/csd.tif"
  ), 
  width = 4.1, height = 6.2, dpi = 600)



sen <- RegionSVGs(sen, grad_cutoff = 0.1, FDR_level = 0.01)


p5f <- TissueFeaturePlot(sen,
                         gene = "CD3e",
                         size = .4,
                         alpha = .6,
                         curve_axes = TRUE) +
  scale_color_gradient(breaks = c(-1, 1, 3),
                       label = c("-1", "1", "3"),
                       low = "#eeeeee",
                       high = "#630202") +
  theme_light() +
  guides(color = guide_colorbar(barwidth = 1,
                                barheight = 6)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -5),
        legend.box.just = "left",
        aspect.ratio = plt_factor); p5f


ggsave(
  plot = p5f,
  filename = paste0(
    path,
    "SENNA/Fig/5_multiome/LN/cd3e.tif"
  ), 
  width = 4.2, height = 6.2, dpi = 600)


p5g <- TissueFeaturePlot(sen,
                         gene = "CD68",
                         size = .4,
                         alpha = .6,
                         curve_axes = TRUE) +
  scale_color_gradient(breaks = c(0, 5, 9),
                      label = c("0", "5", "9"),
                      low = "#eeeeee",
                      high = "#630202") +
  theme_light() +
  guides(color = guide_colorbar(barwidth = 1,
                                barheight = 6)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -5),
        legend.box.just = "left",
        aspect.ratio = plt_factor); p5g

ggsave(
  plot = p5g,
  filename = paste0(
    path,
    "SENNA/Fig/5_multiome/LN/cd68.tif"
  ), 
  width = 4.2, height = 6.2, dpi = 600)


#saveRDS(sen, paste0(path, "dataset/rds/codex_senna.rds"))
