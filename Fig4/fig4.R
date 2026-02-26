library(SENNA)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(cowplot)
library(colorspace)

# path
path <- "/Volumes/KKZR/"
path <- "D:/"

## Islet
isenproc <- function(senna, 
                     knots){
  senna <- TrimmedCurve(senna, knots, type = "islet")
  senna <- GetCurveParam(senna)
  senna <- TissueRegionation(senna)
  return(senna)
}


# Fig. 4a ----

## Data pre-processing
s1 <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus1/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "thymus_1",
  filter.matrix = TRUE,
  image = NULL)

if(min(s1$nCount_Spatial) == 0) s1 <- subset(s1, nCount_Spatial > 0)
s1 <- SCTransform(s1, assay = "Spatial", verbose = FALSE)


s9 <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus9/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "thymus_9",
  filter.matrix = TRUE,
  image = NULL)

if(min(s9$nCount_Spatial) == 0) s9 <- subset(s9, nCount_Spatial > 0)
s9 <- SCTransform(s9, assay = "Spatial", verbose = FALSE)



s10 <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus10/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "thymus_10",
  filter.matrix = TRUE,
  image = NULL)

if(min(s10$nCount_Spatial) == 0) s10 <- subset(s10, nCount_Spatial > 0)
s10 <- SCTransform(s10, assay = "Spatial", verbose = FALSE)



## SENNA; s1
sen1 <- SENNA_Visium(s1,
                     slice_name = "thymus_1",
                     annotation = FALSE)

#AppDat(sen1, 
#       image_path = paste0(path, "dataset/SNUH/Thymus1/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

ks11 <- read.csv(paste0(path, "dataset/knots/thy/s1/k1.csv"),
                 header = TRUE)
ks12 <- read.csv(paste0(path, "dataset/knots/thy/s1/k2.csv"),
                 header = TRUE)
ks13 <- read.csv(paste0(path, "dataset/knots/thy/s1/k3.csv"),
                 header = TRUE)
ks14 <- read.csv(paste0(path, "dataset/knots/thy/s1/k4.csv"),
                 header = TRUE)
ks15 <- read.csv(paste0(path, "dataset/knots/thy/s1/k5.csv"),
                 header = TRUE)
ks16 <- read.csv(paste0(path, "dataset/knots/thy/s1/k6.csv"),
                 header = TRUE)
ks17 <- read.csv(paste0(path, "dataset/knots/thy/s1/k7.csv"),
                 header = TRUE)
ks18 <- read.csv(paste0(path, "dataset/knots/thy/s1/k8.csv"),
                 header = TRUE)
ks19 <- read.csv(paste0(path, "dataset/knots/thy/s1/k9.csv"),
                 header = TRUE)
ks110 <- read.csv(paste0(path, "dataset/knots/thy/s1/k10.csv"),
                  header = TRUE)
ks111 <- read.csv(paste0(path, "dataset/knots/thy/s1/k11.csv"),
                  header = TRUE)
ks112 <- read.csv(paste0(path, "dataset/knots/thy/s1/k12.csv"),
                  header = TRUE)
ks113 <- read.csv(paste0(path, "dataset/knots/thy/s1/k13.csv"),
                  header = TRUE)


for(i in 1:13){
  knots <- get(paste0("ks1", i))
  assign(x = paste0("sen1_", i),
         value = isenproc(sen1, knots))
}

rm("knots", "ks11", "ks12", "ks13",
   "ks14", "ks15", "ks16", "ks17",
   "ks18", "ks19", "ks110", "ks111",
   "ks112", "ks113", "sen1"); gc()

msen1 <- ConvertMultiSENNA(
  list(sen1_1, sen1_2, sen1_3, sen1_4,
       sen1_5, sen1_6, sen1_7, sen1_8,
       sen1_9, sen1_10, sen1_11, sen1_12, 
       sen1_13))
rm("sen1_1", "sen1_2", "sen1_3", "sen1_4",
   "sen1_5", "sen1_6", "sen1_7", "sen1_8",
   "sen1_9", "sen1_10", "sen1_11", "sen1_12", 
   "sen1_13"); gc()

msen1 <- RegionSVGs(msen1,
                    FDR_level = 0.01,
                    grad_cutoff = 0.1,
                    active = FALSE,
                    direction = rep(-1, 
                                    length(msen1@SENNA)))

## SENNA; s9
sen9 <- SENNA_Visium(s9,
                     slice_name = "thymus_9",
                     annotation = FALSE)
#AppDat(sen9, 
#       image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

ks91 <- read.csv(paste0(path, "dataset/knots/thy/s9/k1.csv"),
                 header = TRUE)
ks92 <- read.csv(paste0(path, "dataset/knots/thy/s9/k2.csv"),
                 header = TRUE)
ks93 <- read.csv(paste0(path, "dataset/knots/thy/s9/k3.csv"),
                 header = TRUE)
ks94 <- read.csv(paste0(path, "dataset/knots/thy/s9/k4.csv"),
                 header = TRUE)
ks95 <- read.csv(paste0(path, "dataset/knots/thy/s9/k5.csv"),
                 header = TRUE)
ks96 <- read.csv(paste0(path, "dataset/knots/thy/s9/k6.csv"),
                 header = TRUE)
ks97 <- read.csv(paste0(path, "dataset/knots/thy/s9/k7.csv"),
                 header = TRUE)
ks98 <- read.csv(paste0(path, "dataset/knots/thy/s9/k8.csv"),
                 header = TRUE)
ks99 <- read.csv(paste0(path, "dataset/knots/thy/s9/k9.csv"), 
                 header = TRUE)


for(i in 1:9){
  knots <- get(paste0("ks9", i))
  assign(x = paste0("sen9_", i),
         value = isenproc(sen9, knots))
}

rm("knots", "ks91", "ks92", "ks93",
   "ks94", "ks95", "ks96", "ks97",
   "ks98", "ks99", "sen9"); gc()

msen9 <- ConvertMultiSENNA(
  list(sen9_1, sen9_2, sen9_3, sen9_4,
       sen9_5, sen9_6, sen9_7, sen9_8,
       sen9_9))
rm("sen9_1", "sen9_2", "sen9_3", "sen9_4",
   "sen9_5", "sen9_6", "sen9_7", "sen9_8",
   "sen9_9"); gc()

msen9 <- RegionSVGs(msen9,
                    FDR_level = 0.01,
                    grad_cutoff = 0.1,
                    active = FALSE,
                    direction = rep(-1, length(msen9@SENNA)))


## SENNA; s10
sen10 <- SENNA_Visium(s10,
                      slice_name = "thymus_10",
                      annotation = FALSE)
#AppDat(sen10, 
#       image_path = paste0(path, "dataset/SNUH/Thymus10/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

ks101 <- read.csv(paste0(path, "dataset/knots/thy/s10/k1.csv"),
                  header = TRUE)
ks102 <- read.csv(paste0(path, "dataset/knots/thy/s10/k2.csv"),
                  header = TRUE)
ks103 <- read.csv(paste0(path, "dataset/knots/thy/s10/k3.csv"),
                  header = TRUE)
ks104 <- read.csv(paste0(path, "dataset/knots/thy/s10/k4.csv"),
                  header = TRUE)
ks105 <- read.csv(paste0(path, "dataset/knots/thy/s10/k5.csv"),
                  header = TRUE)
ks106 <- read.csv(paste0(path, "dataset/knots/thy/s10/k6.csv"),
                  header = TRUE)
ks107 <- read.csv(paste0(path, "dataset/knots/thy/s10/k7.csv"),
                  header = TRUE)
ks108 <- read.csv(paste0(path, "dataset/knots/thy/s10/k8.csv"),
                  header = TRUE)
ks109 <- read.csv(paste0(path, "dataset/knots/thy/s10/k9.csv"),
                  header = TRUE)
ks1010 <- read.csv(paste0(path, "dataset/knots/thy/s10/k10.csv"),
                  header = TRUE)

for(i in 1:10){
  knots <- get(paste0("ks10", i))
  assign(x = paste0("sen10_", i),
         value = isenproc(sen10, knots))
}

rm("knots", "ks101", "ks102", "ks103",
   "ks104", "ks105", "ks106", "ks107",
   "ks108", "ks109", "ks1010"); gc()

msen10 <- ConvertMultiSENNA(
  list(sen10_1, sen10_2, sen10_3, sen10_4,
       sen10_5, sen10_6, sen10_7, sen10_8,
       sen10_9, sen10_10))
rm("sen10_1", "sen10_2", "sen10_3", "sen10_4",
   "sen10_5", "sen10_6", "sen10_7", "sen10_8",
   "sen10_9", "sen10_10"); gc()

msen10 <- RegionSVGs(msen10,
                     FDR_level = 0.01,
                     grad_cutoff = 0.1,
                     active = FALSE,
                     direction = rep(-1, length(msen10@SENNA)))


## Fig 4.a; s1
dp <- select(msen1@SENNA[[1]]@Coord[["Spatial"]], region)
dp$region <- dp$region * msen1@SENNA[[2]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[3]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[4]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[5]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[6]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[7]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[8]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[9]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[10]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[11]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[12]]@Coord[["Spatial"]]$region
dp$region <- dp$region * msen1@SENNA[[13]]@Coord[["Spatial"]]$region

dq <- rbind(
  filter(msen1@SENNA[[1]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[2]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[3]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[4]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[5]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[6]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[7]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[8]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[9]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[10]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[11]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[12]]@Coord[["Spatial"]], region == -1),
  filter(msen1@SENNA[[13]]@Coord[["Spatial"]], region == -1)
)

p4a1 <- ggplot() +
  geom_point(aes(X1, X2, color = abs(distance)),
             data = dq,
             alpha = 1,
             size = 1) +
  scale_color_gradient(low = "#bec1b0", 
                       high = "#4b5226",
                       guide = guide_colorbar(
                         ticks = FALSE,
                         breaks = c(min(abs(dq$distance)),
                                    max(abs(dq$distance))),
                         labels = c("Low", "High")
                       )) +
  geom_point(aes(X1, X2), color = "#dddddd",
             data = msen1@SENNA[[1]]@Coord[["Spatial"]][dp$region > 0,],
             alpha = 1,
             size = 1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[1]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[2]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[3]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[4]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[5]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[6]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[7]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[8]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[9]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[10]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[11]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[12]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen1@SENNA[[13]]), 
             color = "#393939", size = .1) +
  theme_test() +
  theme(axis.title = element_blank()) +
  guides(color = "none") + 
  lims(x = c(0, 1), y = c(0, 1)); p4a1

ggsave(plot = p4a1,
       filename = paste0(path, "SENNA/Fig/4_thy/ca1.tif"), 
       height = 3, width = 3, dpi = 600)


## Fig 4.a; s9
dp1 <- select(msen9@SENNA[[1]]@Coord[["Spatial"]], region)
dp1$region <- dp1$region * msen9@SENNA[[2]]@Coord[["Spatial"]]$region
dp1$region <- dp1$region * msen9@SENNA[[3]]@Coord[["Spatial"]]$region
dp1$region <- dp1$region * msen9@SENNA[[4]]@Coord[["Spatial"]]$region
dp1$region <- dp1$region * msen9@SENNA[[5]]@Coord[["Spatial"]]$region
dp1$region <- dp1$region * msen9@SENNA[[6]]@Coord[["Spatial"]]$region
dp1$region <- dp1$region * msen9@SENNA[[7]]@Coord[["Spatial"]]$region
dp1$region <- dp1$region * msen9@SENNA[[8]]@Coord[["Spatial"]]$region
dp1$region <- dp1$region * msen9@SENNA[[9]]@Coord[["Spatial"]]$region

dq1 <- rbind(
  filter(msen9@SENNA[[1]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[2]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[3]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[4]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[5]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[6]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[7]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[8]]@Coord[["Spatial"]], region == -1),
  filter(msen9@SENNA[[9]]@Coord[["Spatial"]], region == -1)
)

p4a2 <- ggplot() +
  geom_point(aes(X1, X2, color = abs(distance)),
             data = dq1,
             alpha = 1,
             size = 1) +
  scale_color_gradient(low = "#bec1b0", 
                       high = "#4b5226",
                       guide = guide_colorbar(
                         ticks = FALSE,
                         breaks = c(min(abs(dq$distance)),
                                    max(abs(dq$distance))),
                         labels = c("Low", "High")
                       )) +
  geom_point(aes(X1, X2), color = "#dddddd",
             data = msen9@SENNA[[1]]@Coord[["Spatial"]][dp1$region > 0,],
             alpha = 1,
             size = 1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[1]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[2]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[3]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[4]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[5]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[6]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[7]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[8]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen9@SENNA[[9]]), 
             color = "#393939", size = .1) +
  theme_test() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank()) +
  guides(color = "none") + 
  lims(x = c(0, 1), y = c(0, 1)); p4a2

ggsave(plot = p4a2,
       filename = paste0(path, "SENNA/Fig/4_thy/ca9.tif"), 
       height = 3, width = 2.73, dpi = 600)



## Fig 4.a; s10
dp2 <- select(msen10@SENNA[[1]]@Coord[["Spatial"]], region)
dp2$region <- dp2$region * msen10@SENNA[[2]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[3]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[4]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[5]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[6]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[7]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[8]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[9]]@Coord[["Spatial"]]$region
dp2$region <- dp2$region * msen10@SENNA[[10]]@Coord[["Spatial"]]$region

dq2 <- rbind(
  filter(msen10@SENNA[[1]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[2]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[3]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[4]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[5]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[6]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[7]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[8]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[9]]@Coord[["Spatial"]], region == -1),
  filter(msen10@SENNA[[10]]@Coord[["Spatial"]], region == -1)
) %>%
  mutate(distance = abs(distance))

p4a3 <- ggplot() +
  geom_point(aes(X1, X2, color = distance),
             data = dq2,
             alpha = 1,
             size = 1) +
  scale_color_gradient(low = "#bec1b0", 
                       high = "#4b5226",
                       name = "Distance",
                       breaks = c(min(dq2$distance + 0.01), 
                                  max(dq2$distance) - 0.01),
                       labels = c("Low", "High")) +
  guides(color = 
           guide_colorbar(
             ticks = FALSE,
             theme = theme(
               legend.ticks = element_blank()))) +
  geom_point(aes(X1, X2), color = "#dddddd",
             data = msen10@SENNA[[1]]@Coord[["Spatial"]][dp2$region > 0,],
             alpha = 1,
             size = 1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[1]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[2]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[3]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[4]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[5]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[6]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[7]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[8]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[9]]), 
             color = "#393939", size = .1) +
  geom_point(aes(X1, X2),
             data = crvtrjry(msen10@SENNA[[10]]), 
             color = "#393939", size = .1) +
  lims(x = c(0, 1), y = c(0, 1)) +
  theme_test() + 
  theme(axis.title = element_blank(),
        axis.text.y = element_blank()); p4a3

ggsave(plot = p4a3,
       filename = paste0(path, "SENNA/Fig/4_thy/ca10.tif"), 
       height = 3, width = 3.47, dpi = 600)



# Fig. 4b ----
msr1 <- msen1@msR
msr9 <- msen9@msR
msr10 <- msen10@msR

set.seed(123)
p4b1 <- RegionVolPlot(msr1, nrepel = 8L, 
                      FDR_level = 0.05,
                      dot_size = 0.8,
                      dot_alpha = 0.8) +
  theme_test() +
  scale_x_continuous(labels = function(x) x^2) +
  ggtitle("S1") +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0))
p4b2 <- RegionVolPlot(msr9, nrepel = 8L, 
                      FDR_level = 0.05,
                      dot_size = 0.8,
                      dot_alpha = 0.8) +
  theme_test() +
  scale_x_continuous(labels = function(x) x^2) +
  ggtitle("S9") +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0))
p4b3 <- RegionVolPlot(msr10, 
                      FDR_level = 0.05,
                      dot_size = 0.8,
                      dot_alpha = 0.8) +
  theme_test() +
  ggtitle("S10") +
  scale_x_continuous(labels = function(x) x^2) +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0))

p4b <- plot_grid(p4b1, p4b2, p4b3,
                 ncol = 1, align = "v")

p4b <- plot_grid(p4b,
                 ggdraw() + 
                   draw_label("Gradient", hjust = 0.4, size = 13),
                 ncol = 1, rel_heights = c(1, 0.04))

p4b <- plot_grid(ggdraw() + 
                   draw_label("Adjusted p-value (-log)", 
                              angle = 90, size = 13),
                 p4b,
                 ncol = 2, rel_widths = c(0.07, 1)); p4b

ggsave(plot = p4b,
       filename = paste0(path, "SENNA/Fig/4_thy/volplot.tif"),
       width = 3.5, height = 3.5*1.7, dpi = 600)


rm(list = ls()); gc()


# Fig. 4c-e ----

# path
path <- "/Volumes/KKZR/"
path <- "D:/"

## Prog
psenproc <- function(senna, 
                     knots){
  senna <- TrimmedCurve(senna, knots, type = "spline")
  senna <- GetCurveParam(senna)
  return(senna)
}


## pre-proc
s1 <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus1/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "thymus_1",
  filter.matrix = TRUE,
  image = NULL)

if(min(s1$nCount_Spatial) == 0) s1 <- subset(s1, nCount_Spatial > 0)
s1 <- SCTransform(s1, assay = "Spatial", verbose = FALSE)

s9 <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus9/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "thymus_9",
  filter.matrix = TRUE,
  image = NULL)

if(min(s9$nCount_Spatial) == 0) s9 <- subset(s9, nCount_Spatial > 0)
s9 <- SCTransform(s9, assay = "Spatial", verbose = FALSE)

# s1
sen1 <- SENNA_Visium(s1,
                     slice_name = "thymus_1",
                     annotation = FALSE)

#AppDat(sen1, 
#       reference_value = "Annotation",
#       image_path = paste0(path, "dataset/SNUH/Thymus1/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

p11 <- read.csv(paste0(path, "dataset/knots/thy/s1/prog/p1.csv"),
                header = TRUE)
p12 <- read.csv(paste0(path, "dataset/knots/thy/s1/prog/p2.csv"),
                header = TRUE)

msen1 <- ConvertMultiSENNA(
  list(psenproc(sen1, p11),
       psenproc(sen1, p12)))

msen1 <- ProgSVGs(msen1,
                  interval = c(0.08, 0.09),
                  FDR_level = 0.01,
                  grad_cutoff = 0.1,
                  active = FALSE)

ProgVolPlot(msen1@msR, 
            FDR_level = 0.01, grad_cutoff = 0.1)


# s9
sen9 <- SENNA_Visium(s9,
                     slice_name = "thymus_9",
                     annotation = TRUE)

#AppDat(sen9, 
#       reference_value = "Annotation",
#       image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

p91 <- read.csv(paste0(path, "dataset/knots/thy/s9/prog/p1.csv"),
                header = TRUE)
p92 <- read.csv(paste0(path, "dataset/knots/thy/s9/prog/p2.csv"),
                header = TRUE)
p93 <- read.csv(paste0(path, "dataset/knots/thy/s9/prog/p3.csv"),
                header = TRUE)

msen9 <- ConvertMultiSENNA(
  list(psenproc(sen9, p91),
       psenproc(sen9, p92),
       psenproc(sen9, p93)))

msen9 <- ProgSVGs(msen9,
                  interval = c(0.055, 0.1, 0.07),
                  FDR_level = 0.01,
                  grad_cutoff = 0.1,
                  active = FALSE)

ProgVolPlot(msen9@msR, 
            FDR_level = 0.01, grad_cutoff = 0.1)


## Fig. 4e
sen11 <- msen1@SENNA[[1]]
sen12 <- msen1@SENNA[[2]]
sen91 <- msen9@SENNA[[1]]
sen92 <- msen9@SENNA[[2]]
sen93 <- msen9@SENNA[[3]]

ps1 <- msen1@msR@Variable_gene[["positive"]]
ns1 <- msen1@msR@Variable_gene[["negative"]]
ps9 <- msen9@msR@Variable_gene[["positive"]]
ns9 <- msen9@msR@Variable_gene[["negative"]]

pi <- intersect(ps1, ps9)
ni <- intersect(ns1, ns9)
rm("ps1", "ps9", "ns1", "ns9"); gc()

genes <- union(pi, ni)

## Pattern Plot
sen11 <- GetCurveParam(sen11)
sen12 <- GetCurveParam(sen12)
sen91 <- GetCurveParam(sen91)
sen92 <- GetCurveParam(sen92)
sen93 <- GetCurveParam(sen93)

sk1 <- Make_simplesktS(sen11, 
                       genelist = genes,
                       nbins = 5,
                       interval = .08)
sk2 <- Make_simplesktS(sen12, 
                       genelist = genes,
                       nbins = 5,
                       interval = .09)
sk3 <- Make_simplesktS(sen91, 
                       genelist = genes,
                       nbins = 5,
                       interval = .055)
sk4 <- Make_simplesktS(sen92, 
                       genelist = genes,
                       nbins = 5,
                       interval = .1)
sk5 <- Make_simplesktS(sen93, 
                       genelist = genes,
                       nbins = 5,
                       interval = .07)

dfheatmap <- function(sskt) {
  df <- data.frame(sskt@sketch)
  df <- df %>%
    mutate(Gene = rownames(df)) %>%
    tidyr::pivot_longer(cols = starts_with("x"),
                        names_to = "CP",
                        values_to = "Expression") %>%
    mutate(
      Group = case_when(
        Gene %in% pi ~ "INC",
        Gene %in% ni ~ "DEC"
      )) %>%
    mutate(Group = factor(Group, 
                          levels = c("INC", "DEC")),
           CP = as.integer(sub("X", "", CP)))
  return(df)
}

df1 <- dfheatmap(sk1) %>% mutate(Dataset = "S1_CA1")
df2 <- dfheatmap(sk2) %>% mutate(Dataset = "S1_CA2")
df3 <- dfheatmap(sk3) %>% mutate(Dataset = "S9_CA1")
df4 <- dfheatmap(sk4) %>% mutate(Dataset = "S9_CA2")
df5 <- dfheatmap(sk5) %>% mutate(Dataset = "S9_CA3")

df <- bind_rows(df1, df2, df3, df4, df5)
rm("df1", "df2", "df3", "df4", "df5"); gc()


p4e <- ggplot(aes(Gene, CP, fill = Expression), 
              data = df) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(divergingx_hcl(100, palette = "RdBu"))) + 
  theme_minimal() +
  facet_grid(cols = vars(Group), 
             rows = vars(Dataset),
             scales = "free", space = "free") +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1), 
        strip.text.x = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = NA,
                                        color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.spacing = unit(0.3, "lines"),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.spacing = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -6)) +
  labs(fill = "Expr") +
  scale_y_reverse(); p4e

p4e <- ggplot_gtable(ggplot_build(p4e))
stript <- which(grepl('strip-t', p4e$layout$name))
facet_colors <- c("#fbaf40", "#a790c0")
k <- 1

for (i in stript) {
  j <- which(grepl('rect', p4e$grobs[[i]]$grobs[[1]]$childrenOrder))
  p4e$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- facet_colors[k]
  k <- k+1
}

grid.draw(p4e)

ggsave(filename = paste0(path, "SENNA/Fig/4_thy/prog_pattern.tif"),
       plot = grid.arrange(p4e),
       width = 10, height = 4.5, dpi = 600)



## Fig. 4c-d
s1ratio <- diff(range(msen1@SENNA[[1]]@Coord$Spatial$X1)) / 
  diff(range(msen1@SENNA[[1]]@Coord$Spatial$X2))

p4c1 <- ShowBins(sen11, bins = 5, interval = .08, 
                 colors = c("#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5"),
                 bg_dot_size = 1.5,
                 dot_size = 1.5) +
  ggtitle("S1_CA1") +
  guides(color = "none") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 1)); p4c1

ggsave(filename = paste0(path, "SENNA/Fig/4_thy/bins11.tif"),
       plot = p4c1,
       width = 3, height = 3 / s1ratio, dpi = 600)

p4c2 <- ShowBins(sen12, bins = 5, interval = .09, 
                 colors = c("#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5"),
                 bg_dot_size = 1.5,
                 dot_size = 1.5) +
  ggtitle("S1_CA2") +
  guides(color = "none") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 1))

ggsave(filename = paste0(path, "SENNA/Fig/4_thy/bins12.tif"),
       plot = p4c2,
       width = 3, height = 3 / s1ratio, dpi = 600)

p4d1 <- ShowBins(sen91, bins = 5, interval = .055, 
                 colors = c("#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5"),
                 bg_dot_size = 1.5,
                 dot_size = 1.5) +
  ggtitle("S9_CA1") +
  guides(color = "none") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 1))

ggsave(filename = paste0(path, "SENNA/Fig/4_thy/bins91.tif"),
       plot = p4d1,
       width = 3, height = 3, dpi = 600)

p4d2 <- ShowBins(sen92, bins = 5, interval = .1, 
                 colors = c("#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5"),
                 bg_dot_size = 1.5,
                 dot_size = 1.5) +
  ggtitle("S9_CA2") +
  guides(color = "none") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 1))

ggsave(filename = paste0(path, "SENNA/Fig/4_thy/bins92.tif"),
       plot = p4d2,
       width = 3, height = 3, dpi = 600)


p4d3 <- ShowBins(sen93, bins = 5, interval = .07, 
                 colors = c("#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5"),
                 bg_dot_size = 1.5,
                 dot_size = 1.5) +
  ggtitle("S9_CA3") +
  labs(color = "CP") +
  guides(color = guide_legend(
    override.aes = list(size = 3)
  )) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 1))

legend <- {
  g <- ggplotGrob(p4d3)
  lid <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[lid]]
}
ggsave(filename = paste0(path, "SENNA/Fig/4_thy/binsleg.tif"),
       plot = legend,
       width = 0.4, height = 1.5, dpi = 600)

p4d3 <- ShowBins(sen93, bins = 5, interval = .07, 
                 colors = c("#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5"),
                 bg_dot_size = 1.5,
                 dot_size = 1.5) +
  ggtitle("S9_CA3") +
  guides(color = "none") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 1)); p4d3

ggsave(filename = paste0(path, "SENNA/Fig/4_thy/bins93.tif"),
       plot = p4d3,
       width = 3, height = 3, dpi = 600)


## EnrichR
library(enrichR)
mdec <- unique(filter(df, Group == "DEC")[["Gene"]])
msen_msig <- enrichr(mdec, databases = "MSigDB_Hallmark_2020")$MSigDB_Hallmark_2020

write.csv(msen_msig, "./Fig/4_thy/msenstab.csv")
