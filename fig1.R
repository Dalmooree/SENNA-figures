library(SENNA) 
library(Seurat)
library(shiny)
library(dplyr)
library(ggplot2)
library(patchwork)



# path
path <- "/Volumes/KKZR/"
path <- "D:/"



#### Fig 1. a-1 ----

surt <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus9/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "thymus_9",
  filter.matrix = TRUE,
  image = NULL)

if(min(surt$nCount_Spatial) == 0) surt <- subset(surt, nCount_Spatial > 0)
surt <- SCTransform(surt, assay = "Spatial", verbose = FALSE)
surt <- RunPCA(surt, assay = "SCT", npcs = 50, verbose = FALSE)
surt <- FindNeighbors(surt, reduction = "pca", 
                      dims = 1:10, verbose = FALSE)
surt <- FindClusters(surt, resolution = 0.8, 
                     algorithm = 4, random.seed = 1, verbose = FALSE)
surt <- RunUMAP(surt, dims = 1:10, verbose = FALSE)
SpatialDimPlot(surt, image.alpha = 0.2)
id <- c("Cortex_intermediate_1",
        "Cortex_posterior_1",
        "Medulla_1",
        "Medulla_3",
        "Cortex_intermediate_2",
        "Cortex_anterior_1",
        "Cortex_posterior_2",
        "Medulla_2",
        "Cortex_anterior_2",
        "Cortex_posterior_3"
)


sen <- SENNA_Visium(surt,
                    slice_name = "thymus_9",
                    annotation = TRUE)



colorset <- c("#7a7a7a", "#7a7a7a", "#fa9b99", "#fa9b99",
              "#7a7a7a", "#7a7a7a", "#7a7a7a", "#fa9b99",
              "#7a7a7a", "#7a7a7a")

#AppDat(sen, colors = colorset, reference_value = "Annotation", image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"))
#knot_picker()

rm(surt); gc()


#### Fig 1.a, prog----
knots <- read.csv(
  paste0(path, "SENNA/Fig/1_main/knots/prog.csv"))

senp <- TrimmedCurve(sen, knots, type = "spline")

pr1 <- ShowCurve(senp, 
                 order_label = FALSE,
                 color_reference = "Annotation",
                 bg_dot_size = 1,
                 bg_dot_alpha = 0.6,
                 knots_size = 1.5,
                 line_size = 1,
                 color = colorset,
                 line_color = "#393939",
                 knots_color = "#000000") + 
  ggrepel::geom_text_repel(aes(X1, X2),
                            data = knots,
                            label = seq(1L, nrow(knots), by = 1L),
                            point.padding = 0.1,
                            box.padding = 0.2,
                            size = 5,
                            seed = 789) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); pr1

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/prog_1.tif"),
       plot = pr1,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")

senp <- GetCurveParam(senp)
pr2 <- CurveParamValid(senp,
                       seed = 111,
                       bg_dot_size = 1,
                       bg_dot_alpha = 0.3,
                       lead_size = 1,
                       lead_type = 2) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); pr2
ggsave(paste0(path,
              "SENNA/Fig/1_main/a/prog_2.tif"),
       plot = pr2,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")


coord <- senp@Coord[["Spatial"]] %>%
  dplyr::mutate(box = ifelse(distance <= 0.1, 1, 0))


pr3 <- ggplot2::ggplot() +
  ggplot2::geom_point(aes(X1, X2, col = t), 
                      size = 1,
                      alpha = .5,
                      data = dplyr::filter(coord,
                                           box == 1)) +
  ggplot2::scale_colour_gradient2(low = "skyblue",
                                  high = "#1b74bc") +
  ggplot2::geom_point(aes(X1, X2),
                      alpha = 0.6, size = 1, col = "lightgray",
                      data = dplyr::filter(coord,
                                           box != 1)) +
  ggplot2::geom_point(aes(X1, X2),
                      alpha = 0.6, size = 1, col = "lightgray",
                      data = dplyr::filter(coord,
                                           is.na(box))) +
  ggplot2::geom_point(aes(X1, X2), data = crvtrjry(senp), 
                      color = "#000000", size = 0.2) +
  ggplot2::geom_point(aes(X1, X2),
                      data = senp@CurveAxis[["knots"]], 
                      color = "#000000", size = 1.2) + 
  ggrepel::geom_text_repel(aes(X1, X2),
                           data = knots,
                           label = seq(1L, nrow(knots), by = 1L),
                           point.padding = 0.1,
                           box.padding = 0.2,
                           size = 5,
                           seed = 789) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); pr3

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/prog_3.tif"),
       plot = pr3,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")

senp <- ProgSVGs(senp,
                 FDR_level = 0.01,
                 active = FALSE,
                 interval = 0.12,
                 grad_cutoff = 0.7)
op1 <- ProgVolPlot(senp,
                   FDR_level = 0.01,
                   dot_size = 0.2,
                   dot_alpha = 0.5) +
  geom_vline(xintercept = sqrt(0.7), lty = "twodash", alpha = 0.7) +
  geom_vline(xintercept = -sqrt(0.7), lty = "twodash", alpha = 0.7) +
  ggplot2::theme_test() +
  ggplot2::labs(y = "p.adj, -log",
                x = "Grad, scaled"); op1
ggsave(paste0(path,
              "SENNA/Fig/1_main/a/op1.tif"),
       plot = op1,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")

op2 <- TissueFeaturePlot(senp,
                         senp@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]][2],
                         size = 1) +
  theme(legend.position = "left",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()) +
  labs(color = "Expr"); op2

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/op2.tif"),
       plot = op2,
       width = 3.45,
       height = 3,
       dpi = 600,
       units = "in")

sskt <- Make_simplesktS(senp,
                        genelist = 
                          c(senp@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]][1:5], senp@Gene[["P.SVGs"]][["Variable_gene"]][["negative"]][1:5]))
rownames(sskt@sketch) <- c("SVG 4", "SVG 2", "SVG 1", "SVG 5", "SVG 3",
                           "SVG 6", "SVG 9", "SVG 8", "SVG 7", "SVG 10")

op3 <- PatternPlot(sskt,
            xlab = FALSE,
            ylab = TRUE,
            legend = FALSE)

svg(paste0(path,
           "SENNA/Fig/1_main/a/op3.svg"),
    width = 6 / 2.54,
    height = 6 / 2.54)
op3
dev.off()


dns <- senp@Coord[["Spatial"]] %>%
  mutate(CLS = senp@Gene$Reference$Annotation)
dns <- dns[complete.cases(dns),]

library(ggridges)
library(scales)
op4 <- ggplot() +
  geom_density_ridges(aes(t, 
                          y = CLS, fill = CLS), 
                      data = dns, alpha = 0.8, color = "#444444",
                      panel_scaling = unit(0.05, "lines")) + 
  scale_fill_manual(values = ggsci::pal_npg(c("nrc"))(10)) + 
  labs(x = "CP",
       y = "Clusters") +
  theme_test() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#dddddd"),
        strip.text = element_text(color = "#000000")); op4

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/op4.tif"),
       plot = op4,
       width = 3,
       height = 3,
       dpi = 600,
       units = "in")


rm("knots", "senp", "coord"); gc()




#### Fig 1.a, islet----

#AppDat(sen, reference_value = "Annotation", colorset)
#knot_picker()

knots <- read.csv(
  paste0(path, "SENNA/Fig/1_main/knots/islet.csv"))

seni <- TrimmedCurve(sen, knots, type = "islet")
is1 <- ShowCurve(seni, 
                 order_label = FALSE,
                 color_reference = "Annotation",
                 bg_dot_alpha = 0.7,
                 bg_dot_size = 1,
                 knots_size = 1.5,
                 line_size = .1,
                 color = colorset,
                 line_color = "#393939",
                 knots_color = "#000000") + 
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); is1

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/isl_1.tif"),
       plot = is1,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")

seni <- GetCurveParam(seni)
seni <- TissueRegionation(seni)

dat <- subset(seni@Coord[["Spatial"]], region == -1)
dat1 <- subset(seni@Coord[["Spatial"]], region == 1)
fun <- seni@CurveAxis[["fun"]]
griddat <- crvtrjry(seni)

set.seed(123123123)
ind <- sample(1:nrow(dat), 6)
d <- lapply(ind, function(i){
  subdat <- dat[i, ] %>%
    dplyr::mutate(index = i)
  
  subdat2 <- data.frame(X1 = trplug_coef(subdat$t, fun$x.coef),
                        X2 = trplug_coef(subdat$t, fun$y.coef),
                        index = i)
  
  subdat <- subdat %>% dplyr::select(X1, X2, index)
  
  sd = BiocGenerics::rbind(subdat, subdat2)
  return(sd)
})



is2 <- ggplot()+
  ggplot2::geom_point(aes(X1, X2), data = dat1, col = 'lightgrey', alpha = .2, size = 1) +
  ggplot2::geom_point(aes(X1, X2), data = dat, col = 'darkgrey', alpha = .2, size = 1) +
  ggplot2::geom_point(aes(X1, X2), data = griddat, color = 'skyblue', size = 0.3) +
  ggplot2::geom_line(aes(X1, X2), linetype = 2,
                     linewidth = 1, col = "#b14743", data = d[[1]]) +
  ggplot2::geom_line(aes(X1, X2), linetype = 2,
                     linewidth = 1, col = "#44729d", data = d[[2]]) +
  ggplot2::geom_line(aes(X1, X2, col = factor(index)), linetype = 2,
                     linewidth = 1, col = "#539045", data = d[[3]]) +
  ggplot2::geom_line(aes(X1, X2, col = factor(index)), linetype = 2,
                     linewidth = 1, col = "#7e5d55", data = d[[4]]) +
  ggplot2::geom_line(aes(X1, X2, col = factor(index)), linetype = 2,
                     linewidth = 1, col = "#8e73ae", data = d[[5]]) +
  ggplot2::geom_line(aes(X1, X2, col = factor(index)), linetype = 2,
                     linewidth = 1, col = "#d48640", data = d[[6]]) +
  ggplot2::geom_point(aes(X1, X2),
                      shape = 21, size = 1 * 1.2,
                      col = "#b14743", fill = 'grey',
                      data = d[[1]]) +
  ggplot2::geom_point(aes(X1, X2),
                      shape = 21, size = 1 * 1.2,
                      col = "#44729d", fill = 'grey',
                      data = d[[2]]) +
  ggplot2::geom_point(aes(X1, X2),
                      shape = 21, size = 1 * 1.2,
                      col = "#539045", fill = 'grey',
                      data = d[[3]]) +
  ggplot2::geom_point(aes(X1, X2),
                      shape = 21, size = 1 * 1.2,
                      col = "#7e5d55", fill = 'grey',
                      data = d[[4]]) +
  ggplot2::geom_point(aes(X1, X2),
                      shape = 21, size = 1 * 1.2,
                      col = "#8e73ae", fill = 'grey',
                      data = d[[5]]) +
  ggplot2::geom_point(aes(X1, X2),
                      shape = 21, size = 1 * 1.2,
                      col = "#d48640", fill = 'grey',
                      data = d[[6]]) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); is2

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/isl_2.tif"),
       plot = is2,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")
rm("dat", "dat1", "d", "fun", "griddat", "ind"); gc()


is3 <- ggplot2::ggplot() +
  ggplot2::geom_point(aes(X1, X2, col = distance), 
                      alpha = 0.7, size = 1,
                      data = dplyr::filter(seni@Coord[["Spatial"]],
                                           distance <= 0)) +
  ggplot2::scale_colour_gradient(high = "#dec9b3", 
                                 low = "#ca6804") +
  ggplot2::geom_point(aes(X1, X2),
                      alpha = 0.6, size = 1, col = "lightgray",
                      ,
                      data = dplyr::filter(seni@Coord[["Spatial"]],
                                           distance > 0)) +
  ggplot2::geom_point(aes(X1, X2), data = crvtrjry(seni), 
                      color = "black", size = 0.1) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); is3

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/isl_3.tif"),
       plot = is3,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")

rm(seni); gc()




#### Fig 1.a, regio----

#AppDat(sen, reference_value = "Annotation", colorset)
#knot_picker()

knots <- read.csv(
  paste0(path, "SENNA/Fig/1_main/knots/regio.csv"))

ranges <- ggplot_build(pr1)$layout$panel_params[[1]]

senr <- FullCurve(sen, knots, type = "spline")
rg1 <- ShowCurve(senr, 
                 order_label = FALSE,
                 color_reference = "Annotation",
                 bg_dot_alpha = 0.7,
                 bg_dot_size = 1,
                 knots_size = 1.5,
                 line_size = .1,
                 color = colorset,
                 line_color = "#393939",
                 knots_color = "#000000") +
  coord_cartesian(xlim = ranges$x.range, 
                  ylim = ranges$y.range,
                  expand = FALSE) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); rg1

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/rg_1.tif"),
       plot = rg1,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")

senr <- GetCurveParam(senr)
senr <- TissueRegionation(senr)

rg2 <- CurveParamValid(senr,
                       seed = 321,
                       bg_dot_size = 1,
                       bg_dot_alpha = 0.3,
                       lead_size = 1,
                       lead_type = 2) +
  coord_cartesian(xlim = ranges$x.range, 
                  ylim = ranges$y.range,
                  expand = FALSE) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", 
                                       color = "black", linewidth = 1)); rg2

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/rg_2.tif"),
       plot = rg2,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")

rg3 <- ShowCSDistance(senr,
                      dot_size = 1,
                      dot_alpha = .7,
                      line_size = .2,
                      high = "#8e73ae",
                      medium = "#eeeeee",
                      low = "#32a02c") +
  coord_cartesian(xlim = ranges$x.range, 
                  ylim = ranges$y.range,
                  expand = FALSE) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1)); rg3

ggsave(paste0(path,
              "SENNA/Fig/1_main/a/rg_3.tif"),
       plot = rg3,
       width = 6,
       height = 6,
       dpi = 600,
       units = "cm")




#### Fig 1.b, mSENNA ----
surt <- Load10X_Spatial(data.dir = paste0(path, 
                                          "dataset/SNUH/Thymus9/outs"),
                        filename = "filtered_feature_bc_matrix.h5",
                        assay = "Spatial",
                        slice = "thymus_9",
                        filter.matrix = TRUE,
                        image = NULL)

if(min(surt$nCount_Spatial) == 0) surt <- subset(surt, nCount_Spatial > 0)
surt <- SCTransform(surt, assay = "Spatial", verbose = FALSE)
surt <- RunPCA(surt, assay = "SCT", npcs = 50, verbose = FALSE)
surt <- FindNeighbors(surt, reduction = "pca", 
                      dims = 1:10, verbose = FALSE)
surt <- FindClusters(surt, resolution = 0.8, 
                     algorithm = 4, random.seed = 1, verbose = FALSE)
surt <- RunUMAP(surt, dims = 1:10, verbose = FALSE)
SpatialDimPlot(surt, image.alpha = 0.2)
id <- c("Cortex_intermediate_1",
        "Cortex_posterior_1",
        "Medulla_1",
        "Medulla_3",
        "Cortex_intermediate_2",
        "Cortex_anterior_1",
        "Cortex_posterior_2",
        "Medulla_2",
        "Cortex_anterior_2",
        "Cortex_posterior_3"
)


sen <- SENNA_Visium(surt,
                    slice_name = "thymus_9",
                    annotation = TRUE)



colorset_b <- c("#7a7a7a", "#7a7a7a", "#fa9b99", "#fa9b99",
                "#7a7a7a", "#7a7a7a", "#7a7a7a", "#fa9b99",
                "#7a7a7a", "#7a7a7a")


#AppDat(sen, reference_value = "Annotation", colorset_b)
#knot_picker()



knots1 <- read.csv(paste0(path, "dataset/knots/thy/s9/m1.csv"), 
                   header = T)
knots2 <- read.csv(paste0(path, "dataset/knots/thy/s9/m2.csv"), 
                   header = T)
knots3 <- read.csv(paste0(path, "dataset/knots/thy/s9/m3.csv"), 
                   header = T)
knots4 <- read.csv(paste0(path, "dataset/knots/thy/s9/m4.csv"), 
                   header = T)
knots5 <- read.csv(paste0(path, "dataset/knots/thy/s9/m5.csv"), 
                   header = T)
knots6 <- read.csv(paste0(path, "dataset/knots/thy/s9/m6.csv"), 
                   header = T)
knots7 <- read.csv(paste0(path, "dataset/knots/thy/s9/m7.csv"), 
                   header = T)
knots8 <- read.csv(paste0(path, "dataset/knots/thy/s9/m8.csv"), 
                   header = T)
knots9 <- read.csv(paste0(path, "dataset/knots/thy/s9/m9.csv"), 
                   header = T)


sen1 <- TrimmedCurve(sen, knots1, type = "islet")
sen2 <- TrimmedCurve(sen, knots2, type = "islet")
sen3 <- TrimmedCurve(sen, knots3, type = "islet")
sen4 <- TrimmedCurve(sen, knots4, type = "islet")
sen5 <- TrimmedCurve(sen, knots5, type = "islet")
sen6 <- TrimmedCurve(sen, knots6, type = "islet")
sen7 <- TrimmedCurve(sen, knots7, type = "islet")
sen8 <- TrimmedCurve(sen, knots8, type = "islet")
sen9 <- TrimmedCurve(sen, knots9, type = "islet")


griddat1 <- crvtrjry(sen1)
griddat2 <- crvtrjry(sen2)
griddat3 <- crvtrjry(sen3)
griddat4 <- crvtrjry(sen4)
griddat5 <- crvtrjry(sen5)
griddat6 <- crvtrjry(sen6)
griddat7 <- crvtrjry(sen7)
griddat8 <- crvtrjry(sen8)
griddat9 <- crvtrjry(sen9)
val <- sen1@Coord[["Spatial"]]
val <- merge(x = val,
             y = data.frame(sen1@Gene[["Reference"]][["Annotation"]]),
             by = "row.names",
             all = TRUE)
rownames(val) <- val[, 1]
val <- val[, 2:4]
colnames(val) <- c("X1", "X2", "CLS")

b2 <- ggplot2::ggplot() +
  ggplot2::geom_point(aes(X1, X2, color = CLS),
                      data = val,
                      alpha = 0.6,
                      size = .8) +
  ggplot2::scale_color_manual(values = colorset_b) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat1, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat2, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat3, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat4, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat5, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat6, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat7, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat8, color = "#393939", size = .1) +
  ggplot2::geom_point(aes(X1, X2),
                      data = griddat9, color = "#393939", size = .1) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", 
                                       color = "black", linewidth = 1)); b2

ggsave(paste0(path,
              "SENNA/Fig/1_main/b/_1.tif"),
       plot = b2,
       width = 6,
       height = 6,
       device = "tif",
       units = "cm",
       dpi = 600)


rm("griddat1", "griddat2", "griddat3",
   "griddat4", "griddat5", "griddat6",
   "griddat7", "griddat8", "griddat9",
   "knots1", "knots2", "knots3",
   "knots4", "knots5", "knots6",
   "knots7", "knots8", "knots9", 
   "id", "val", "surt"); gc()



b31 <- ShowCurve(sen1, order_label = FALSE, cluster_label = TRUE,
                 bg_dot_size = .7, bg_dot_alpha = .7,
                 knots_size = 0, knots_color = "#393939",
                 line_color = "#393939",
                 color = colorset_b) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1))

ggsave("./Figure/1/b/_31.svg",
       plot = b31,
       width = 6,
       height = 6,
       device = "svg",
       units = "cm")

b32 <- ShowCurve(sen2, order_label = FALSE, cluster_label = TRUE,
                 bg_dot_size = .7, bg_dot_alpha = .7,
                 knots_size = 0, knots_color = "#393939",
                 line_color = "#393939",
                 color = colorset_b) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1))

ggsave("./Figure/1/b/_32.svg",
       plot = b32,
       width = 6,
       height = 6,
       device = "svg",
       units = "cm")

b33 <- ShowCurve(sen3, order_label = FALSE, cluster_label = TRUE,
                 bg_dot_size = .7, bg_dot_alpha = .7,
                 knots_size = .1, knots_color = "#393939",
                 line_color = "#393939",
                 color = colorset_b) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1))

ggsave("./Figure/1/b/_33.svg",
       plot = b33,
       width = 6,
       height = 6,
       device = "svg",
       units = "cm")

b34 <- ShowCurve(sen4, order_label = FALSE, cluster_label = TRUE,
                 bg_dot_size = .7, bg_dot_alpha = .7,
                 knots_size = 0, knots_color = "#393939",
                 line_color = "#393939",
                 color = colorset_b) +
  ggplot2::theme_void() +
  theme(legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.background = element_rect(fill = "white", color = "black", linewidth = 1))

ggsave("./Figure/1/b/_34.svg",
       plot = b34,
       width = 6,
       height = 6,
       device = "svg",
       units = "cm")



sen1 <- GetCurveParam(sen1); sen1 <- TissueRegionation(sen1)
sen1 <- ActiveIdent(sen1, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen2 <- GetCurveParam(sen2); sen2 <- TissueRegionation(sen2)
sen2 <- ActiveIdent(sen2, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen3 <- GetCurveParam(sen3); sen3 <- TissueRegionation(sen3)
sen3 <- ActiveIdent(sen3, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen4 <- GetCurveParam(sen4); sen4 <- TissueRegionation(sen4)
sen4 <- ActiveIdent(sen4, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen5 <- GetCurveParam(sen5); sen5 <- TissueRegionation(sen5)
sen5 <- ActiveIdent(sen5, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen6 <- GetCurveParam(sen6); sen6 <- TissueRegionation(sen6)
sen6 <- ActiveIdent(sen6, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen7 <- GetCurveParam(sen7); sen7 <- TissueRegionation(sen7)
sen7 <- ActiveIdent(sen7, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen8 <- GetCurveParam(sen8); sen8 <- TissueRegionation(sen8)
sen8 <- ActiveIdent(sen8, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))
sen9 <- GetCurveParam(sen9); sen9 <- TissueRegionation(sen9)
sen9 <- ActiveIdent(sen9, c("Medulla_1", "Medulla_2", "Medulla_3", "Medulla_4"))

msen <- ConvertMultiSENNA(list(sen1, sen2, sen3,
                               sen4, sen5, sen6,
                               sen7, sen8, sen9))
rm("sen1", "sen2", "sen3",
   "sen4", "sen5", "sen6",
   "sen7", "sen8", "sen9"); gc()

msen <- RegionSVGs(msen,
                   FDR_level = 0.05,
                   grad_cutoff = 0.01,
                   direction = rep(-1, 9),
                   active = TRUE)
#saveRDS(msen, "./Figure/1/b/msenna.rds")
b4 <- RegionVolPlot(msen,
                    FDR_level = 0.05,
                    dot_size = .1,
                    dot_alpha = .5) +
  ggplot2::xlim(-2.5, 2.5) +
  ggplot2::theme_test() +
  ggplot2::theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank()) +
  ggplot2::labs(y = "p.adj, -log",
                x = "Grad, scaled"); b4 

ggsave("./Figure/1/b/_4.svg",
       plot = b4,
       width = 5.5,
       height = 6,
       device = "svg",
       units = "cm")
