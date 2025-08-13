library(SENNA)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(scales)
library(ggridges)

# path
path <- "/Volumes/KKZR/"
path <- "D:/"
#

# S1. h-polynomial ----

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

sen9 <- SENNA_Visium(s9,
                     slice_name = "thymus_9",
                     annotation = FALSE)
#AppDat(sen9, 
#       image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

k1 <- read.csv(paste0(path, "SENNA/Supplementary/knots/k1.csv"),
               header = TRUE)


sen9 <- FullCurve(sen9, k1, type = "spline")

k1 <- mutate(k1, N = 1:4)
set.seed(3);id <- sample(1:nrow(sen9@Coord[["Spatial"]]), 1)

ps1 <- ShowCurve(sen9, 
                 order_label = FALSE,
                 bg_dot_size = 1.5,
                 bg_dot_alpha = .3) +
  ggrepel::geom_text_repel(aes(X1, X2, label = paste0("t=", N)), 
                           data = k1,
                           seed = 1) + 
  geom_point(aes(X1, X2), data = sen9@Coord[["Spatial"]][id,],
             size = 1.5, color = "darkred", shape = 1, fill = NA) +
  labs(x = "X", y = "Y") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))); ps1

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S1_hpol.png"),
       plot = ps1,
       width = 6, height = 6, dpi = 1000)




# S3. LUAD density, full ----
 
senr1 <- readRDS(paste0(path, "dataset/rds/xp_sen_r1.rds"))
senr2 <- readRDS(paste0(path, "dataset/rds/xp_sen_r2.rds"))

R1 <- mutate(senr1@Coord[["Spatial"]],
             type = senr1@Gene[["Reference"]][["Annotation"]])
R2 <- mutate(senr2@Coord[["Spatial"]],
             type = senr2@Gene[["Reference"]][["Annotation"]])

ct <- unique(as.character(c(R1[["type"]], R2[["type"]])))

preprocess <- function(df, source_label) {
  df %>%
    filter(type %in% ct) %>%
    mutate(
      distance = rescale(abs(distance), to = c(0, 1)),
      Source = source_label)}

R1 <- preprocess(R1, "R1")
R2 <- preprocess(R2, "R2")

ridg <- bind_rows(R1, R2)


ps3 <- 
  ggplot(ridg, aes(x = distance, y = type, fill = type)) +
  geom_density_ridges(alpha = 0.8, scale = 0.7) +
  facet_wrap(~ Source, ncol = 2) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Distance from curve axis (scaled)") +
  theme_light() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#dddddd"),
        strip.text = element_text(color = "#000000")); ps3

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S3_luad.png"),
       plot = ps3,
       width = 9, height = 10, dpi = 1000)

rm(list = setdiff(ls(), "path")); gc()

# S Table 1 ~ 3 (gene list) ----

## LUAD VHD
senv1 <- readRDS(paste0(path, "dataset/rds/VHD8_r1.rds"))
senv2 <- readRDS(paste0(path, "dataset/rds/VHD8_r2.rds"))

v1 <- senv1@Gene[["R.SVGs"]][["Report"]]
v1 <- v1 %>%
  filter(Gene %in%
           c(senv1@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
           senv1@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])) %>%
  mutate(Sample = "R1")

v2 <- senv2@Gene[["R.SVGs"]][["Report"]]
v2 <- v2 %>%
  filter(Gene %in%
           c(senv2@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
             senv2@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])) %>%
  mutate(Sample = "R2")

vhd <- rbind(v1, v2) %>%
  mutate(Sample = as.factor(Sample))
write.csv(vhd, paste0(path, "SENNA/Supplementary/stable/stab1_vhd.csv"))

rm("senv1", "senv2", "v1", "v2", "vhd") ; gc()




## Thymus islet

isenproc <- function(senna, 
                     knots){
  senna <- TrimmedCurve(senna, knots, type = "islet")
  senna <- GetCurveParam(senna)
  senna <- TissueRegionation(senna)
  return(senna)
}

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


sen1 <- SENNA_Visium(s1,
                     slice_name = "thymus_1",
                     annotation = FALSE)

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

sen9 <- SENNA_Visium(s9,
                     slice_name = "thymus_9",
                     annotation = FALSE)

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

sen10 <- SENNA_Visium(s10,
                      slice_name = "thymus_10",
                      annotation = FALSE)
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
   "ks108", "ks109", "ks1010", "sen10"); gc()

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

mi1 <- msen1@msR@Report
mi1 <- mi1 %>%
  filter(Gene %in%
           c(msen1@msR@Variable_gene[["positive"]],
             msen1@msR@Variable_gene[["negative"]])) %>%
  mutate(Sample = "S1")

mi9 <- msen9@msR@Report
mi9 <- mi9 %>%
  filter(Gene %in%
           c(msen9@msR@Variable_gene[["positive"]],
             msen9@msR@Variable_gene[["negative"]])) %>%
  mutate(Sample = "S9")

mi10 <- msen10@msR@Report
mi10 <- mi10 %>%
  filter(Gene %in%
           c(msen10@msR@Variable_gene[["positive"]],
             msen10@msR@Variable_gene[["negative"]])) %>%
  mutate(Sample = "S10")

mseni <- rbind(mi1, mi9, mi10) %>%
  mutate(Sample = as.factor(Sample))
write.csv(mseni, 
          paste0(path, "SENNA/Supplementary/stable/stab2_thy_islet.csv"))

rm(list = setdiff(ls(), c("s1", "s9", "path"))); gc()




## Thymus progression
psenproc <- function(senna, 
                     knots){
  senna <- TrimmedCurve(senna, knots, type = "spline")
  senna <- GetCurveParam(senna)
  return(senna)
}

sen1 <- SENNA_Visium(s1,
                     slice_name = "thymus_1",
                     annotation = FALSE)

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

sen9 <- SENNA_Visium(s9,
                     slice_name = "thymus_9",
                     annotation = TRUE)

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

cgp <- intersect(msen1@msR@Variable_gene[["positive"]],
                 msen9@msR@Variable_gene[["positive"]])
cgn <- intersect(msen1@msR@Variable_gene[["negative"]],
                 msen9@msR@Variable_gene[["negative"]])
mp1 <- msen1@msR@Report
mp1 <- mp1 %>%
  filter(Gene %in% c(cgp, cgn)) %>%
  mutate(Sample = "S1")

mp9 <- msen9@msR@Report
mp9 <- mp9 %>%
  filter(Gene %in% c(cgp, cgn)) %>%
  mutate(Sample = "S9")

msenp <- rbind(mp1, mp9) %>%
  mutate(Sample = as.factor(Sample))
write.csv(msenp, 
          paste0(path, "SENNA/Supplementary/stable/stab3_thy_prog.csv"))

rm(list = ls()); gc()




## CODEX regionation
path <- "/Volumes/KKZR/"
path <- "D:/"

sen <- readRDS(paste0(path, "dataset/rds/codex_senna.rds"))

c1 <- sen@Gene[["R.SVGs"]][["Report"]]
c1 <- c1 %>%
  filter(Gene %in%
           c(sen@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
             sen@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]]))

write.csv(c1, 
          paste0(path, "SENNA/Supplementary/stable/stab4_codex.csv"))




# S2. Piecewise linear curve axis ----

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

sen9 <- SENNA_Visium(s9,
                     slice_name = "thymus_9",
                     annotation = FALSE)
#AppDat(sen9, 
#       image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

l1 <- read.csv(paste0(path, "SENNA/Supplementary/knots/lin1.csv"),
               header = TRUE)


ls9 <- FullCurve(sen9, l1, type = "spline")
ll9 <- FullCurve(sen9, l1, type = "straight")

ps1 <- ShowCurve(ls9, 
                 order_label = FALSE,
                 bg_dot_size = 1,
                 bg_dot_alpha = .5,
                 line_size = 1) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()); ps1

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S2_splca.tiff"),
       plot = ps1,
       width = 4, height = 4, dpi = 600)

ps2 <- ShowCurve(ll9, 
                 order_label = FALSE,
                 bg_dot_size = 1,
                 bg_dot_alpha = .5)+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()); ps2

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S2_linca.tiff"),
       plot = ps2,
       width = 4, height = 4, dpi = 600)

ls9 <- GetCurveParam(ls9)
ls9 <- TissueRegionation(ls9)
ll9 <- GetCurveParam(ll9)
ll9 <- TissueRegionation(ll9)

ps3 <- ShowRegions(ls9,
                   dot_size = 1,
                   dot_alpha = .5); ps3

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S2_regspl.tiff"),
       plot = ps3,
       width = 4.6, height = 4, dpi = 300)

ps4 <- ShowRegions(ll9,
                   dot_size = 1,
                   dot_alpha = .5); ps4

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S2_reglin.tiff"),
       plot = ps4,
       width = 4.6, height = 4, dpi = 300)

ps5 <- ShowRegions(ls9,
                   dot_size = 13,
                   dot_alpha = .5,
                   line_size = .6) +
  guides(color = "none") + 
  coord_cartesian(xlim = c(0.25, 0.30),
                  ylim = c(0.25, 0.30)); ps5

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S2_regspl_loop.tiff"),
       plot = ps5,
       width = 4, height = 4, dpi = 300)


ps6 <- ShowRegions(ll9,
                   dot_size = 13,
                   dot_alpha = .5,
                   line_size = .6) +
  guides(color = "none") + 
  coord_cartesian(xlim = c(0.25, 0.30),
                  ylim = c(0.25, 0.30)); ps6

ggsave(filename = paste0(path, "SENNA/Supplementary/sfig/S2_reglin_loop.tiff"),
       plot = ps6,
       width = 4, height = 4, dpi = 300)



