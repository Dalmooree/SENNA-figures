## SENNA; Fig. 3
  
library(SENNA)
library(Seurat)
library(spacexr)
library(enrichR)
library(ggplot2)
library(dplyr)
library(ggridges)
library(scales)
library(ggsci)
library(arrow)

# set path
path <- getwd()

##

# Xen5k; pre-proc ----

object <- LoadXenium(data.dir = 
                       paste0(path, 
                              "dataset/PublicXenium/Xenium_Prime_Human_Lung_Cancer_FFPE_outs/"),
                     assay = "Xenium",
                     fov = "fov")
object <- subset(object, subset = nCount_Xenium > 0)

object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object,
                               verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)

object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch",
  verbose = FALSE
  )

gc()

DefaultAssay(object) <- "sketch"
object <- ScaleData(object,
                    verbose = FALSE)
object <- RunPCA(object, 
                 assay = "sketch", 
                 reduction.name = "pca.sketch",
                 verbose = FALSE)
object <- FindNeighbors(object, 
                        reduction = "pca.sketch", 
                        dims = 1:50,
                        verbose = FALSE)
object <- RunUMAP(object, 
                  reduction = "pca.sketch", 
                  reduction.name = "umap.sketch", 
                  return.model = T, 
                  dims = 1:50,
                  verbose = FALSE)
gc()

## RCTD
ref <- readRDS(
  paste0(path,
         "dataset/Reference/droplet_normal_lung_blood_filtered/vhd1/luad_ref_P2_d1.rds"))

cluster <- as.factor(ref$free_annotation)
counts <- ref[["RNA"]]$counts
cluster <- droplevels(cluster)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
rctd_ref <- spacexr::Reference(counts, cluster, nUMI, min_UMI = 1)
cntxe <- object[["Xenium"]]$counts
cells <- colnames(object[["Xenium"]])
coord <- GetTissueCoordinates(object)
rownames(coord) <- coord$cell
coord <- coord[,1:2]
query <- spacexr::SpatialRNA(coords = coord, 
                             counts = cntxe, 
                             nUMI = colSums(cntxe))

rm("cntxe", "nUMI", "cluster", "counts", "cells", "coord"); gc()

rctd <- create.RCTD(query, 
                    rctd_ref, 
                    max_cores = parallel::detectCores()-1)
rctd <- run.RCTD(rctd,
                 doublet_mode = "doublet")
object <- AddMetaData(object, metadata = rctd@results$results_df)
rm("rctd_ref", "ref", "query") ; gc()
rm("Q_mat", "K_val", "N_X", "X_vals", "S_mat", "SQ_mat"); gc()

object$first_type <- as.character(object$first_type)
object$first_type[is.na(object$first_type)] <- "Unknown"
object <- ProjectData(
  object = object,
  assay = "Xenium",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type"),
  verbose = FALSE)

gc()


## Build Niches
DefaultAssay(object) <- "Xenium"
object <- BuildNicheAssay(object = object, 
                          fov = "fov",
                          group.by = "full_first_type",
                          niches.k = 5, 
                          neighbors.k = 30)

ImageDimPlot(object, 
             fov = "fov",
             group.by = "niches", 
             size = 0.3, 
             dark.background = F) +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))

saveRDS(object, paste0(path, "dataset/rds/xp_p2d1_nch.rds"))

## Fig. 3b ----

object <- readRDS(paste0(path, "dataset/rds/xp_p2d1_nch.rds"))
DefaultAssay(object) <- "Xenium"

xenrotate <- function(coords) {
  x <- coords[,1]; y <- coords[,2]
  coords[,1] <- max(y) - y
  coords[,2] <- x - min(x)
  return(coords)
  }

object@images$fov@boundaries$centroids@coords <- 
  xenrotate(object@images$fov@boundaries$centroids@coords)

egfr <- FetchData(object, vars = "EGFR")

p3b <- ImageFeaturePlot(object, "EGFR",
                        alpha = 0.3,
                        size = 0.2,
                        dark.background = FALSE,
                        min.cutoff = min(egfr),
                        max.cutoff = max(egfr),
                        axes = TRUE,
                        coord.fixed = FALSE) +
  scale_fill_gradientn(colors = viridis::turbo(n = 10),
                       values = rescale(c(
                         min(egfr),
                         median(unique(egfr[[1]])) - sd(egfr[[1]]),
                         median(unique(egfr[[1]])),
                         median(unique(egfr[[1]])) + sd(egfr[[1]]),
                         max(egfr)
                         )),
                       breaks = c(round(min(egfr), 0), 
                                  round(median(unique(egfr[[1]])), 0),
                                  6)) + 
  theme_light() +
  labs(color = "EGFR") +
  guides(fill = guide_colorbar(barwidth = .8,
                               barheight = 3)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        legend.ticks = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, -2, 0, -3)); p3b

ggsave(plot = p3b,
       filename = paste0(path, "SENNA/Fig/3_luad/ncs_new/egfr.tif"),
       width = 6.14/2, height = 3.35/2, dpi = 600)

# Xen5k; SENNA ----
object <- readRDS(paste0(path, "dataset/rds/xp_p2d1_nch.rds"))
DefaultAssay(object) <- "Xenium"
Idents(object) <- object$full_first_type
sen <- SENNA_Xenium(object, 
                    fov = "fov",
                    annotation = TRUE)
sen <- AddReference(sen, "Niches", as.factor(object$niches))

rm("object", "xenrotate", "egfr"); gc()


## Pick knots
#AppDat(sen, reference_value = "Niches", colorset = c("#8e73ae", "#44729d", "#d48640", "#539045", "#b14743"))
#knot_picker()


## CA
xenk1 <- read.csv(
  paste0(path, "dataset/knots/xp_knots_r1.csv"), 
  header = TRUE)
xenk2 <- read.csv(
  paste0(path, "dataset/knots/xp_knots_r2.csv"), 
  header = TRUE)


## Fig. 3c ----
senr1 <- TrimmedCurve(senna = sen, knot_df = xenk1, type = "islet")
senr2 <- TrimmedCurve(senna = sen, knot_df = xenk2, type = "islet")

p3c <- ShowCurve(senr1, 
                 color_reference = "Niches",
                 colors = c("#b14743", "#a790c0", "#1b74bc", "#67b665","#fbaf40"),
                 line_color = "#000000", 
                 knots_color = "#000000",
                 bg_dot_size = 0.1, 
                 bg_dot_alpha = 0.1,
                 knots_size = 0,
                 line_size = .8,
                 order_label = FALSE) + 
  geom_path(aes(X1, X2), 
             data = crvtrjry(senr2), 
             color = "#000000", 
             size = .8) +
  theme_light() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, .7),
                  expand = TRUE) +
  labs(color = "Niches") +
  guides(color = 
           guide_legend(override.aes = list(size = 1.5,
                                            alpha = 1),
                        title.hjust = 1.5)) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6,
                                     margin = margin(l=2, unit = "pt")),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -5),
          legend.spacing.y = unit(-2, units = "pt"),
          legend.box.just = "left",
          legend.key.height = unit(10, "pt")); p3c
  
  
ggsave(plot = p3c,
       filename = paste0(path, "SENNA/Fig/3_luad/ncs_new/xenca.tif"),
       width = 6.14/2, height = 3.35 / 2 , dpi = 600)
  
  
# Cell type density ----
senr1 <- GetCurveParam(senr1)
senr1 <- TissueRegionation(senr1)
saveRDS(senr1, paste0(path, "dataset/rds/xp_sen_r1.rds"))
rm(senr1); gc()

senr2 <- GetCurveParam(senr2)
senr2 <- TissueRegionation(senr2)
saveRDS(senr2, paste0(path, "dataset/rds/xp_sen_r2.rds"))
rm(senr2); gc()

  
senr1 <- readRDS(paste0(path, "dataset/rds/xp_sen_r1.rds"))
senr2 <- readRDS(paste0(path, "dataset/rds/xp_sen_r2.rds"))
  
R1 <- mutate(senr1@Coord[["Spatial"]],
             type = senr1@Gene[["Reference"]][["Annotation"]])
R2 <- mutate(senr2@Coord[["Spatial"]],
             type = senr2@Gene[["Reference"]][["Annotation"]])
  
ct <- c(
  "CD8+ Memory-Effector T",
  "CD4+ Memory-Effector T",
  "Myeloid Dendritic Type 1",
  "Macrophage",
  "Natural Killer",
  "Capillary",
  "Pericyte",
  "Alveolar Epithelial Type 2"
  )
  
  
preprocess <- function(df, source_label) {
  df <- df %>%
    filter(type %in% ct) %>%
    mutate(
      distance = rescale(abs(distance), to = c(0, 1)),
      Source = source_label)
  df$type <- factor(df$type, levels = rev(ct))
  return(df)}
  
R1 <- preprocess(R1, "R1")
R2 <- preprocess(R2, "R2")
  
ridg <- bind_rows(R1, R2)
  
  
## Fig. 3d ----
p3d <- 
  ggplot(ridg, aes(x = distance, y = type, fill = type)) +
  geom_density_ridges(alpha = 0.6, scale = .7, linewidth = .3) +
  scale_fill_manual(values = pal_npg("nrc")(10)) + 
  facet_wrap(~ Source, ncol = 2) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Distance from curve axis (scdaled)") +
  theme_light() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        strip.background = element_rect(fill = "#dfdfdf"),
        strip.text = element_text(color = "#000000", size = 7,
                                  margin = margin(t=2, b=2))) ; p3d

rm("senr1", "senr2", "R1", "R2", "ct"); gc()

ggsave(
  plot = p3d,
  filename = paste0(path, "SENNA/Fig/3_luad/ncs_new/xen_density1.tif"), 
  width = 9.3/2, height = 6/2,  dpi = 600)
  
  
  
# VHD; pre-proc ----
object <- Load10X_Spatial(data.dir = paste0(path,
                                            "dataset/PublicVHD/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment1_spatial/"),
                          bin.size = 8)

object <- subset(object, subset = nCount_Spatial.008um > 3)
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object,
                               verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)

## SENNA
DefaultAssay(object) <- "Spatial.008um"
gc()
sen <- SENNA_Visium(object, 
                    slice_name = "slice1.008um",
                    annotation = FALSE)

sen <- AddReference(sen, 
                    var_name = "EGFR",
                    reference = sen@Gene[["Spatial"]][["EGFR"]])

saveRDS(sen, paste0(path, "dataset/rds/vhd_expr1_008.rds"))
rm(object); gc()


### Pick knots based on EGFR expression level
#AppDat(sen, 
#       reference_value = "EGFR",
#       colorset = c("#dddddd", "red"))
#knot_picker()

## R1
vhdk1 <- read.csv(
  paste0(path,
         "dataset/knots/vhd8_knots_1.csv"), 
  header = TRUE)

senv1 <- TrimmedCurve(senna = sen, knot_df = vhdk1, type = "islet")
senv1 <- GetCurveParam(senv1)
senv1 <- TissueRegionation(senv1)
ShowRegions(senv1)
ShowCSDistance(senv1, direction = -1)
senv1 <- RegionSVGs(senv1,
                    FDR_level = 0.01, 
                    grad_cutoff = 0.1,
                    direction = -1)

saveRDS(senv1, 
        paste0(path, "dataset/rds/VHD8_r1.rds"))

## R2
vhdk2 <- read.csv(
  paste0(path,
         "dataset/knots/vhd8_knots_2.csv"), 
  header = TRUE)

senv2 <- TrimmedCurve(senna = sen, knot_df = vhdk2, type = "islet")
senv2 <- GetCurveParam(senv2)
senv2 <- TissueRegionation(senv2)
ShowRegions(senv2)
senv2 <- RegionSVGs(senv2,
                    FDR_level = 0.01, 
                    grad_cutoff = 0.1,
                    direction = -1)

saveRDS(senv2, 
        paste0(path, "dataset/rds/VHD8_r2.rds"))

## Fig.3e ----
senv1 <- readRDS(paste0(path, "dataset/rds/VHD8_r1.rds"))
senv2 <- readRDS(paste0(path, "dataset/rds/VHD8_r2.rds"))



## Fig 3e
mxd <- abs(min(c(senv1@Coord[["Spatial"]][["distance"]],
                 senv2@Coord[["Spatial"]][["distance"]])))
bks <- c(mxd * 0.2,
         mxd * 0.8)

p3e <- ggplot() +
  geom_point(aes(X1, X2, col = abs(distance)), 
             alpha = 0.8, size = 0.1,
             data = dplyr::filter(senv1@Coord[["Spatial"]],
                                  distance <= 0)) +
  geom_point(aes(X1, X2, col = abs(distance)), 
             alpha = 0.8, size = 0.1,
             data = dplyr::filter(senv2@Coord[["Spatial"]],
                                  distance <= 0)) +
  scale_colour_gradient(low = "#dec9b3", 
                        high = "#ca6804",
                        name = "Distance",
                        limits = c(0, mxd),
                        breaks = bks,
                        labels = c("low", "high"),
                        expand = c(0, 0)) +
  geom_point(aes(X1, X2),
             alpha = 0.5, size = 0.1, col = "#dddddd",
             data = intersect(
               filter(senv1@Coord[["Spatial"]],
                      distance > 0)[,c("X1", "X2")],
               filter(senv2@Coord[["Spatial"]],
                      distance > 0)[,c("X1", "X2")]
               )) +
  geom_path(aes(X1, X2), data = crvtrjry(senv1), 
             color = "#000000", size = .8) +
  geom_path(aes(X1, X2), data = crvtrjry(senv2), 
             color = "#000000", size = .8) +
  theme_light() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.margin = margin(0, -2, 0, -3),
        legend.ticks = element_blank()) +
  guides(color = guide_colorbar(barwidth = 1,
                                barheight = 4)); p3e

ggsave(
  plot = p3e,
  filename = paste0(path, "SENNA/Fig/3_luad/ncs_new/vhdca1.tif"), 
  width = 7.5 / 2, height = 6.7 / 2, dpi = 600)


# VHD EnrichR (MSigDB_Hallmark_2020) ----
senv1 <- readRDS(paste0(path, "dataset/rds/VHD8_r1.rds"))
senv2 <- readRDS(paste0(path, "dataset/rds/VHD8_r2.rds"))

r1p <- senv1@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]]
r1n <- senv1@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]]
r2p <- senv2@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]]
r2n <- senv2@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]]
p2d <- setdiff(r2p, r1p)
p1d <- setdiff(r1p, r2p)
n2d <- setdiff(r2n, r1n)
n1d <- setdiff(r1n, r2n)

pi <- intersect(r1p, r2p)
ni <- intersect(r1n, r2n)


# Stab1
stab1 <- data.frame(
  "Category" = c(
    "R1-increased (R1I)", 
    "R1-decreased (R1D)", 
    "R2-increased (R2I)", 
    "R2-decreased (R2D)",
    "R1 only-increased (R1I - R2I)",
    "R1 only-decreased (R1D - R2D)",
    "R2 only-increased (R2I - R1I)",
    "R2 only-decreased (R2D - R1D)"
  ),
  "Count" = c(
    length(r1p),
    length(r1n),
    length(r2p),
    length(r2n),
    length(p1d),
    length(n1d),
    length(p2d),
    length(n2d) 
  ),
  stringsAsFactors = FALSE
)

write.csv(stab1,
          file.path(path, "SENNA", "__NCS__", "__Supple___", "SupplementaryTable1_.csv"))


## R1 only (pos)
enpr1f <- enrichr(p1d, databases = "MSigDB_Hallmark_2020")
enpr1 <- enpr1f$MSigDB_Hallmark_2020
ts1 <- filter(enpr1, Adjusted.P.value <= 0.05)[["Term"]]

## R2 only (pos)
enpr2f <- enrichr(p2d, databases = "MSigDB_Hallmark_2020")
enpr2 <- enpr2f$MSigDB_Hallmark_2020
ts2 <- filter(enpr2, Adjusted.P.value <= 0.05)[["Term"]]

## R1 only (neg)
ennr1f <- enrichr(n1d, databases = "MSigDB_Hallmark_2020")
ennr1 <- ennr1f$MSigDB_Hallmark_2020
ts3 <- filter(ennr1, Adjusted.P.value <= 0.05)[["Term"]]

## R2 only (neg)
ennr2f <- enrichr(n2d, databases = "MSigDB_Hallmark_2020")
ennr2 <- ennr2f$MSigDB_Hallmark_2020
ts4 <- filter(ennr2, Adjusted.P.value <= 0.05)[["Term"]]

ts <- union(ts1, union(ts2, union(ts3, ts4)))

enpr1 <- enpr1%>%
  filter(.$Term %in% ts) %>%
  mutate(Source = "R1 only (increased)",
         adjp = Adjusted.P.value) %>%
  select(Term, adjp, Source)

enpr2 <- enpr2 %>%
  filter(.$Term %in% ts) %>%
  mutate(Source = "R2 only (increased)",
         adjp = Adjusted.P.value) %>%
  select(Term, adjp, Source)

ennr1 <- ennr1 %>%
  filter(.$Term %in% ts) %>%
  mutate(Source = "R1 only (decreased)",
         adjp = Adjusted.P.value) %>%
  select(Term, adjp, Source)

ennr2 <- ennr2 %>%
  filter(.$Term %in% ts) %>%
  mutate(Source = "R2 only (decreased)",
         adjp = Adjusted.P.value) %>%
  select(Term, adjp, Source)


## Fig. 3f ----
### Supplementary figure. S2
dat <- rbind(enpr1, enpr2, ennr1, ennr2) %>%
  mutate(Source = as.factor(Source),
         adjp = -log10(adjp))

p3f_full <- ggplot() +
  geom_point(aes(x = Term, 
                 y = adjp,
                 shape = Source, 
                 color = Source),
             size = 2,
             stroke = 1.2,
             fill = NA,
             data = dat) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_light() +
  scale_shape_manual(
    values = c(
      "R1 only (increased)" = 24,
      "R2 only (increased)" = 24,
      "R1 only (decreased)" = 25,
      "R2 only (decreased)" = 25)
    ) +
  scale_color_manual(
    values = c(
      "R1 only (increased)" =  "#8e73ae",
      "R2 only (increased)" = "#539045",
      "R1 only (decreased)" = "#7e5d55" ,
      "R2 only (decreased)" = "#d48640")
    ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) +
  labs(y = "Adjusted p-value (-log10)"); p3f_full

ggsave(plot = p3f_full,
       filename = paste0(path,
                         "SENNA/Supplementary/sfig/S4_enrich.tif"),
       height = 8, width = 10, dpi = 600)


### Fig. 3f
filts <- c("Interferon Gamma Response",
           "Inflammatory Response", 
           "Angiogenesis",
           "Interferon Alpha Response",
           "IL-2/STAT5 Signaling",
           "IL-6/JAK/STAT3 Signaling",
           "Complement")

dafa <- dat %>% filter(Term %in% filts)

p3f <- ggplot() +
  geom_point(aes(x = Term, 
                 y = adjp,
                 shape = Source, 
                 color = Source),
             size = 1.8,
             alpha = 0.8,
             stroke = 1.2,
             fill = NA,
             data = dafa) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_light() +
  scale_shape_manual(
    values = c(
      "R1 only (increased)" = 24,
      "R2 only (increased)" = 24,
      "R1 only (decreased)" = 25,
      "R2 only (decreased)" = 25)
    ) +
  scale_color_manual(
    values = c(
      "R1 only (increased)" =  "#8e73ae",
      "R2 only (increased)" = "#539045",
      "R1 only (decreased)" = "#7e5d55" ,
      "R2 only (decreased)" = "#d48640")
    ) +
  guides(color = guide_legend(keywidth = unit(.5, "lines"),
                              keyheight = unit(.5, "lines"),
                              override.aes = list(size = 1.2,
                                                  stroke = .8)),
         shape = guide_legend(label.position = "right",
                              label.hjust = 0,
                              label.theme = element_text(size = 6, 
                                                         margin = margin(l = -3)))) +
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1,
                                   size = 6),
        axis.text.y =  element_text(size = 6),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title.x = element_blank(),
        legend.margin = margin(t=0, r=0, b=0, l=0),
        legend.box.margin = margin(0, 0, 0, -10)) +
  labs(y = "Adjusted p-value (-log10)"); p3f

ggsave(plot = p3f,
       filename = paste0(path,
                         "SENNA/Fig/3_luad/ncs_new/enrichr1.png"),
       height = 6/2, width = 7.7/2, dpi = 600)


