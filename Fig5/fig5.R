library(SENNA)
library(Seurat)
library(ggplot2)
library(dplyr)

# path
path <- "/Volumes/One Touch/"
path <- "D:/"


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



saveRDS(cod, paste0(path, "dataset/rds/codex_seurat.rds"))


# SENNA
cod <- readRDS(paste0(path, "dataset/rds/codex_seurat.rds"))

sen <- SENNA_CODEX(cod,
                   assay = "Akoya",
                   fov_name = "HBM754.WKLP.262",
                   annotation = TRUE)

TissueFOV(sen, 
          dot_size = .3, 
          dot_alpha = .3, 
          colors = DiscretePalette(
            length(unique(cod$seurat_clusters)), 
            palette = "parade", 
            shuffle = FALSE))

colorset <- DiscretePalette(
  length(unique(cod$seurat_clusters)), 
  palette = "parade", 
  shuffle = FALSE)

#AppDat(sen, reference_value = "Annotation", colorset = colorset)
#knot_picker()

ck <- read.csv(
  paste0(path, "dataset/knots/codex_knots.csv"), 
  header = TRUE)

sen <- FullCurve(sen, ck, "spline")

plt_factor <- (max(sen@Coord$Spatial$X2) - min(sen@Coord$Spatial$X2)) / 
  (max(sen@Coord$Spatial$X1) - min(sen@Coord$Spatial$X1))

p5a <- ShowCurve(sen,
                  color_reference = "Annotation",
                  colors = colorset,
                  line_color = "#000000",
                  knots_color = "#000000",
                  bg_dot_size = .1,
                  bg_dot_alpha = .3,
                  knots_size = .3,
                  line_size = .3,
                  order_label = FALSE) +
  theme_light() +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = .8))) +
  labs(color = "Cluster") +
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5a

get_legend <- function(p) {
  g <- ggplotGrob(p)
  legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[legend_index]]
}

legend <- get_legend(p5a)

ggsave(plot = cowplot::ggdraw(legend),
       filename = paste0(
         path,
         "SENNA/Fig/5_codex/clslegend.tif"
       ), 
       width = .5, height = 3, dpi = 600)

p5a <- ShowCurve(sen,
                 color_reference = "Annotation",
                 colors = colorset,
                 line_color = "#000000",
                 knots_color = "#000000",
                 bg_dot_size = .1,
                 bg_dot_alpha = .3,
                 knots_size = .3,
                 line_size = .3,
                 order_label = FALSE) +
  theme_light() +
  guides(color = "none") +
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5a

ggsave(
  plot = p5a,
  filename = paste0(
    path,
    "SENNA/Fig/5_codex/ca.tif"
  ), 
  height = 6 * plt_factor * .6, width = 6 * .6, dpi = 600)


sen <- GetCurveParam(sen)
sen <- TissueRegionation(sen)

p5b <- ShowCSDistance(sen,
                      dot_size = .1,
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
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5b

legend <- get_legend(p5b)

ggsave(plot = cowplot::ggdraw(legend),
       filename = paste0(
         path,
         "SENNA/Fig/5_codex/csdlegend.tif"
       ), 
       width = .5, height = 2, dpi = 600)


p5b <- ShowCSDistance(sen,
                      dot_size = .1,
                      dot_alpha = .3) +
  scale_color_gradient2(high = "#32a02c",
                        mid = "#dedede",
                        low = "#6a3d9a",
                        breaks = c(min(sen@Coord$Spatial$distance) * .8, 
                                   0, 
                                   max(sen@Coord$Spatial$distance) * .8),
                        labels = c("-", 0, "+")) +
  guides(color = "none") + 
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5b
ggsave(
  plot = p5b,
  filename = paste0(
    path,
    "SENNA/Fig/5_codex/csd.tif"
  ), 
  height = 6 * plt_factor * .6, width = 6 * .6, dpi = 600)



sen <- RegionSVGs(sen, grad_cutoff = 0.1, FDR_level = 0.01)

incg <- sen@Gene$R.SVGs$Variable_gene$positive
decg <- sen@Gene$R.SVGs$Variable_gene$negative

id <- 1
for(g in incg){
  TissueFeaturePlot(sen,
                    gene = g,
                    size = .1,
                    alpha = .4,
                    low = "#eeeeee",
                    curve_axes = TRUE) +
    theme_light() +
    theme(axis.title = element_blank(),
          axis.text = element_blank())
  
  ggsave(paste0(path,
                "SENNA/Fig/5_codex/svgs/INC",
                id, "_", g, ".tif"),
         height = 6 * plt_factor * .6, width = 7 * .6, dpi = 600)
  
  id <- id + 1
}


id <- 1
for(g in decg){
  TissueFeaturePlot(sen,
                    gene = g,
                    size = .1,
                    alpha = .4,
                    low = "#eeeeee",
                    curve_axes = TRUE) +
    theme_light() +
    theme(axis.title = element_blank(),
          axis.text = element_blank())
  
  ggsave(paste0(path,
                "SENNA/Fig/5_codex/svgs/DEC",
                id, "_", g, ".tif"),
         height = 6 * plt_factor * .6, width = 7 * .6, dpi = 600)
  
  id <- id + 1
}


p5c <- TissueFeaturePlot(sen,
                         gene = "CD3e",
                         size = .4,
                         alpha = .3,
                         low = "#eeeeee",
                         curve_axes = TRUE) +
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5c

legend <- get_legend(p5c)

ggsave(plot = cowplot::ggdraw(legend),
       filename = paste0(
         path,
         "SENNA/Fig/5_codex/cd3elegend.tif"
       ), 
       width = .5, height = 1.8, dpi = 600)

p5c <- TissueFeaturePlot(sen,
                         gene = "CD3e",
                         size = .4,
                         alpha = .3,
                         low = "#eeeeee",
                         curve_axes = TRUE) +
  guides(color = "none") +
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5c

ggsave(
  plot = p5c,
  filename = paste0(
    path,
    "SENNA/Fig/5_codex/cd3e.tif"
  ), 
  height = 6 * plt_factor * .6, width = 6 * .6, dpi = 600)


p5d <- TissueFeaturePlot(sen,
                         gene = "CD68",
                         size = .4,
                         alpha = .3,
                         low = "#eeeeee",
                         curve_axes = TRUE) +
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5d

legend <- get_legend(p5d)

ggsave(plot = cowplot::ggdraw(legend),
       filename = paste0(
         path,
         "SENNA/Fig/5_codex/cd68legend.tif"
       ), 
       width = .6, height = 1.8, dpi = 600)

p5d <- TissueFeaturePlot(sen,
                         gene = "CD68",
                         size = .4,
                         alpha = .3,
                         low = "#eeeeee",
                         curve_axes = TRUE) +
  guides(color = "none") +
  theme_light() +
  theme(axis.title = element_blank(),
        axis.text = element_blank()); p5d

ggsave(
  plot = p5d,
  filename = paste0(
    path,
    "SENNA/Fig/5_codex/cd68.tif"
  ), 
  height = 6 * plt_factor * .6, width = 6 * .6, dpi = 600)


saveRDS(sen, paste0(path, "dataset/rds/codex_senna.rds"))

sen <- readRDS(paste0(path, "dataset/rds/codex_senna.rds"))
