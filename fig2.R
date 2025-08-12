library(SENNA)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(parallel)
library(SPARK)
library(Giotto)

# path
path <- "/Volumes/KKZR/"
path <- "D:/"




source(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/ssource.R"))

#buffer


#### Shortest distance----

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
s1 <- RunPCA(s1, assay = "SCT", npcs = 50, verbose = FALSE)
s1 <- FindNeighbors(s1, reduction = "pca", 
                    dims = 1:50, verbose = FALSE)
s1 <- FindClusters(s1, resolution = 2, 
                   algorithm = 4, verbose = FALSE)
s1 <- RunUMAP(s1, dims = 1:50, verbose = FALSE)



## 1. Spline curve

sen <- SENNA_Visium(s1,
                    slice_name = "thymus_1",
                    annotation = TRUE)

#AppDat(sen, 
#       reference_value = "Annotation",
#       image_path = paste0(path, "dataset/SNUH/Thymus1/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

prsim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/pr.csv"),
                  header = TRUE)
clsim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/cl.csv"),
                  header = TRUE)

senc <- TrimmedCurve(sen, clsim, type = "islet")
sen <- FullCurve(sen, prsim, type = "spline")


qry <- sen@Coord$Spatial[rownames(sen@Coord$Spatial) == "TCTAGTTATCAGAAGA-1",]


### Validation in progression analysis
sen <- GetCurveParam(sen)
ShowCurve(sen)

sen <- rescalecp(sen)
bks <- c(0.1, 0.25, 0.6, 1)

bks <- sen@Coord$Spatial[sapply(bks, 
                                function(q) which.min(abs(sen@Coord$Spatial$tprime - q))), ]

bks <- {
  X1 <- X2 <- c()
  id <- 1
  for(t in bks$t){
    X1[id] <- plug_coef(t, sen@CurveAxis$fun$x.coef)
    X2[id] <- plug_coef(t, sen@CurveAxis$fun$y.coef)
    id <- id + 1
  }
  
  data.frame(X1 = X1, X2 = X2, t = bks$t, distance = bks$distance, tprime = bks$tprime, row.names = rownames(bks))
}

qry <- sen@Coord$Spatial[rownames(sen@Coord$Spatial) == rownames(qry),]


ShowCurve(sen, 
          order_label = FALSE, 
          colors = "#dddddd",
          bg_dot_size = 1,
          line_color = ggsci::pal_npg("nrc")(1),
          knots_color = ggsci::pal_npg("nrc")(1),
          knots_size = 0.1) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5,
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[1, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5, 
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[2, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5, 
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[3, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5, 
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[4, c("X1", "X2")])) +
  geom_point(aes(X1, X2), data = bks,
             size = 1.3, col = ggsci::pal_npg("nrc")(1)) +
  geom_point(aes(X1, X2), data = qry, 
             shape = 21, col = ggsci::pal_npg("nrc")(3)[3], size = 1.5,
             fill = ggsci::pal_npg("nrc")(3)[3]) +
  geom_point(aes(X1, X2), 
             data = data.frame(X1 = plug_coef(qry$t, sen@CurveAxis$fun$x.coef), X2 = plug_coef(qry$t, sen@CurveAxis$fun$y.coef)),
             shape = 21, col = ggsci::pal_npg("nrc")(3)[3], size = 1.5,
             fill = ggsci::pal_npg("nrc")(3)[3]) +
  geom_line(aes(X1, X2), linetype = "solid", 
            col = ggsci::pal_npg("nrc")(3)[3], 
            linewidth = 1, 
            data = rbind(qry[,c("X1","X2")], data.frame(X1 = plug_coef(qry$t, sen@CurveAxis$fun$x.coef), X2 = plug_coef(qry$t, sen@CurveAxis$fun$y.coef)))) +
  ggrepel::geom_label_repel(aes(X1, X2, 
                               label = paste0("t=", round(t, 1))), 
                           data = bks,
                           box.padding = 0.2,
                           size = 3.5,
                           fill = alpha("#ffffff", 0.5),
                           color = "#000000",
                           seed = 123) + 
  ggrepel::geom_label_repel(aes(X1, X2, 
                               label = paste0("t=", round(t, 1))), 
                           data = data.frame(X1 = plug_coef(qry$t, sen@CurveAxis$fun$x.coef), 
                                             X2 = plug_coef(qry$t, sen@CurveAxis$fun$y.coef),
                                             t = sen@Coord$Spatial$t[rownames(sen@Coord$Spatial) %in% rownames(qry)]),
                           fontface = "bold",
                           box.padding = 0.5,
                           nudge_y = 0.02,
                           fill = alpha("#ffffff", 0.8),
                           color = "#000000",,
                           label.padding = unit(0.4, "lines"),
                           size = 4.5,
                           seed = 1) +
  labs(x = "X", y = "Y")

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/prog_ca2.tif"),
       dpi = 600, width = 4, height = 0.7*4)


td <- seq(min(sen@Coord$Spatial$t), 
          max(sen@Coord$Spatial$t), length.out = 100)
tx <- sapply(td, FUN = plug_coef, coef_mat = sen@CurveAxis$fun$x.coef)
ty <- sapply(td, FUN = plug_coef, coef_mat = sen@CurveAxis$fun$y.coef)
td <- tibble(t = td, X1 = tx, X2 = ty)
td$Distance <- sqrt((td$X1 - qry$X1)^2 + (td$X2 - qry$X2)^2)

bks <- bks %>%
  mutate(distance = purrr::map_dbl(t, ~ {
    nearest_idx <- which.min(abs(td$t - .x))
    td$Distance[nearest_idx]
  }))


ggplot() +
  geom_vline(xintercept = bks$t[1],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[2],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[3],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[4],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_line(aes(t, Distance), data = td,
            col = ggsci::pal_npg("nrc")(4)[4], linewidth = 1) +
  geom_vline(xintercept = qry$t,
             col = ggsci::pal_npg("nrc")(3)[3],
             linewidth = .8, linetype = "longdash") +
  labs(x = "t") +
  scale_x_continuous(breaks = c(0, bks$t[4]),
                     labels = c("0",
                                format(bks$t[4], digits = 2))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 9)) +
  ggrepel::geom_label_repel(
    aes(x = t, y = distance, label = paste0("t = ", round(t, 1))),
    data = bks,
    box.padding = 0.2,
    size = 3.5,
    fill = alpha("#ffffff", 0.8),
    color = "#000000",
    direction = "y",
    seed = 111
  ) +
  geom_label(
    aes(x = t, y = .5, label = paste0("t = ", round(t, 1))),
    data = qry,
    fontface = "bold",
    fill = alpha("#ffffff", 0.8),
    color = "#000000",,
    label.padding = unit(0.4, "lines"),
    size = 4.5
  )

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/prog_dist2.tif"),
       dpi = 600, width = 4, height = 0.7*4)


### Validation in islet analysis

senc <- GetCurveParam(senc)
ShowCurve(senc)

bks <- seq(from = min(senc@Coord$Spatial$t), 
           to = 10, length.out = 4)

bks <- senc@Coord$Spatial[sapply(bks, 
                                function(q) which.min(abs(senc@Coord$Spatial$t - q))), ]

bks <- {
  X1 <- X2 <- c()
  id <- 1
  for(t in bks$t){
    X1[id] <- trplug_coef(t, senc@CurveAxis$fun$x.coef)
    X2[id] <- trplug_coef(t, senc@CurveAxis$fun$y.coef)
    id <- id + 1
  }
  
  data.frame(X1 = X1, X2 = X2, t = bks$t, row.names = rownames(bks))
}

qry <- senc@Coord$Spatial[rownames(senc@Coord$Spatial) == rownames(qry),]

ShowCurve(senc, 
          order_label = FALSE, 
          colors = "#dddddd",
          bg_dot_size = 1,
          line_color = ggsci::pal_npg("nrc")(1),
          knots_color = ggsci::pal_npg("nrc")(1),
          knots_size = 0.1) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[1, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[2, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[3, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .5, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[4, c("X1", "X2")])) +
  geom_point(aes(X1, X2), data = bks,
             size = 1.3, col = ggsci::pal_npg("nrc")(1)) +
  geom_point(aes(X1, X2), data = qry, 
             shape = 21, fill = ggsci::pal_npg("nrc")(3)[3], 
             col = ggsci::pal_npg("nrc")(3)[3], size = 1.5) +
  geom_point(aes(X1, X2), 
             data = data.frame(X1 = trplug_coef(qry$t, senc@CurveAxis$fun$x.coef), X2 = trplug_coef(qry$t, senc@CurveAxis$fun$y.coef)),
             shape = 21, fill = ggsci::pal_npg("nrc")(3)[3], size = 1.5,
             col = ggsci::pal_npg("nrc")(3)[3]) +
  geom_line(aes(X1, X2), linetype = "solid", 
            col = ggsci::pal_npg("nrc")(3)[3], 
            linewidth = 1, 
            data = rbind(qry[,c("X1","X2")], 
                         data.frame(X1 = trplug_coef(qry$t, senc@CurveAxis$fun$x.coef), 
                                    X2 = trplug_coef(qry$t, senc@CurveAxis$fun$y.coef)))) +
  ggrepel::geom_label_repel(aes(X1, X2, 
                               label = paste0("t=", round(t, 1))), 
                           data = bks,,
                           point.padding = 0.01,
                           box.padding = 0.2,
                           size = 3,
                           fill = alpha("#ffffff", 0.5),
                           color = "#000000",
                           seed = 123) + 
  ggrepel::geom_label_repel(aes(X1, X2, 
                               label = paste0("t=", round(t, 2))), 
                           data = data.frame(X1 = trplug_coef(qry$t, senc@CurveAxis$fun$x.coef), 
                                             X2 = trplug_coef(qry$t, senc@CurveAxis$fun$y.coef),
                                             t = senc@Coord$Spatial$t[rownames(senc@Coord$Spatial) %in% rownames(qry)]),
                           fontface = "bold",
                           point.padding = 0.01,
                           box.padding = 0.1,
                           nudge_y = -0.01,
                           fill = alpha("#ffffff", 0.8),
                           color = "#000000",,
                           label.padding = unit(0.4, "lines"),
                           size = 4.5,
                           seed = 1) +
  labs(x = "X", y = "Y")

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/isl_ca1.tif"),
       dpi = 600, width = 4, height = 0.7*4)

td <- seq(min(senc@Coord$Spatial$t), 
          max(senc@Coord$Spatial$t), length.out = 100)
tx <- sapply(td, FUN = trplug_coef, coef_mat = senc@CurveAxis$fun$x.coef)
ty <- sapply(td, FUN = trplug_coef, coef_mat = senc@CurveAxis$fun$y.coef)
td <- tibble(t = td, X1 = tx, X2 = ty)
td$Distance <- sqrt((td$X1 - qry$X1)^2 + (td$X2 - qry$X2)^2)

bks <- bks %>%
  mutate(distance = purrr::map_dbl(t, ~ {
    nearest_idx <- which.min(abs(td$t - .x))
    td$Distance[nearest_idx]
  }))

ggplot() +
  geom_vline(xintercept = bks$t[1],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[2],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[3],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[4],
             col = "#444444", alpha = 0.5,
             linewidth = .3, linetype = "longdash") +
  geom_line(aes(t, Distance), data = td,
            col = ggsci::pal_npg("nrc")(4)[4], linewidth = 1) +
  geom_vline(xintercept = qry$t,
             col = ggsci::pal_npg("nrc")(3)[3],
             linewidth = .8, linetype = "longdash")+ 
  labs(x = "t") +
  scale_x_continuous(breaks = c(1, max(td)),
                     labels = c("1", 
                                format(max(td), digits = 1))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 9)) +
  ggrepel::geom_label_repel(
    aes(x = t, y = distance, label = paste0("t = ", round(t, 1))),
    data = bks,
    box.padding = 0.2,
    size = 3.5,
    fill = alpha("#ffffff", 0.8),
    color = "#000000",
    direction = "y",
    seed = 123
  ) +
  geom_label(
    aes(x = t, y = .23, label = paste0("t = ", round(t, 2))),
    data = qry,
    fontface = "bold",
    fill = alpha("#ffffff", 0.8),
    color = "#000000",,
    label.padding = unit(0.4, "lines"),
    size = 4.5
  )

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/isl_dist2.tif"),
       dpi = 600, width = 4, height = 0.7*4)


rm(list = ls());gc()

# buffer





#### SVGs detection----

# path
path <- "/Volumes/KKZR/"
path <- "D:/"




source(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/ssource.R"))

s1 <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus9/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Visium",
  filter.matrix = TRUE,
  image = NULL)

if(min(s1$nCount_Spatial) == 0) s1 <- subset(s1, nCount_Spatial > 0)
surt <- s1
s1 <- SCTransform(s1, assay = "Spatial", variable.features.n = 2000, verbose = FALSE)

sen <- SENNA_Visium(s1,
                    slice_name = "Visium")
rm(s1); gc()

AppDat(sen, 
       image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"),
       image_resolution = "lowres")
#knot_picker()



prsim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/svgp.csv"),
                  header = TRUE)
rgsim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/svgr.csv"),
                  header = TRUE)
issim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/svgi.csv"),
                  header = TRUE)

seni <- TrimmedCurve(sen, issim, type = "islet")
senr <- FullCurve(sen, rgsim, type = "spline")
senp <- FullCurve(sen, prsim, type = "spline")

seni <- GetCurveParam(seni)
senr <- GetCurveParam(senr)
senp <- GetCurveParam(senp)

seni <- TissueRegionation(seni)
senr <- TissueRegionation(senr)

pypath <- ifelse(Sys.info()["sysname"] == "Windows",
                 "C:/Users/lmseo/AppData/Local/miniconda3/envs/py11/python.exe",
                 "/opt/miniconda3/envs/py310/bin/python3.10")  

instrs <- createGiottoInstructions(save_plot = FALSE,
                                   show_plot = FALSE,
                                   python_path = pypath)
gio <- createGiottoVisiumObject(
  visium_dir = paste0(path, "dataset/SNUH/Thymus1/outs"),
  expr_data = 'filter',
  png_name = 'tissue_lowres_image.png',
  gene_column_index = 2,
  instructions = instrs)


spkcoord <- sen@Coord[["Spatial"]]
#spkcoord <- read.csv(paste0(path, "/SENNA/Fig/2_sim/simulation_benchmarking/coord.csv"), row.names = 1, check.names = FALSE)



set.seed(1)
baseline <- rnorm(n = 2000, mean = 10, sd = 1)



### prog
library(parallel)
t0 <- Sys.time()
prog <- mclapply(seq(from = 0.1, to = 1, length.out = 30),
                 progsimulation,
                 MoreArgs = list(
                   n = 1,
                   poipar = baseline,
                   knot = prsim,
                   seurat = surt,
                   senna = senp,
                   giotto = gio,
                   sparkarg = spkcoord,
                   pval = 0.05),
                 mc.preschedule = TRUE,
                 mc.cores = 1)
prog <- do.call(rbind, prog)
t0 <- Sys.time() - t0
#saveRDS(prog, "./Figure/benchmark/sim_report/prog.rds")
#write.csv(prog, "./Figure/benchmark/sim_report/prog.csv")


# ver Windows
t0 <- Sys.time()
prog <- lapply(seq(from = 0.1, to = 1, length.out = 30),
               progsimulation,
               n = 1,
               poipar = baseline,
               knot = prsim,
               seurat = surt,
               senna = senp,
               giotto = gio,
               sparkarg = spkcoord,
               pval = 0.05)
prog <- do.call(rbind, prog)
t0 <- Sys.time() - t0
saveRDS(prog, "./Figure/2_sim/benchmark/prog.rds")
write.csv(prog, "./Figure/2_sim/benchmark/prog.csv")



### regio

t1 <- Sys.time()
regio <- mclapply(seq(from = 0.1, to = 1, length.out = 30),
                  regiosimulation,
                  MoreArgs = list(
                    n = 1,
                    poipar = baseline,
                    knot = rgsim,
                    seurat = surt,
                    senna = senr,
                    giotto = gio,
                    sparkarg = spkcoord,
                    pval = 0.05),
                  mc.preschedule = TRUE,
                  mc.cores = 5)
regio <- do.call(rbind, regio)
t1 <- Sys.time() - t1
#saveRDS(regio, "./Figure/benchmark/sim_report/regio.rds")
#write.csv(regio, "./Figure/benchmark/sim_report/regio.csv")


# ver Windows
t1 <- Sys.time()
regio <- lapply(seq(from = 0.1, to = 1, length.out = 30),
                regiosimulation,
                n = 1,
                poipar = baseline,
                knot = rgsim,
                seurat = surt,
                senna = senr,
                giotto = gio,
                sparkarg = spkcoord,
                pval = 0.05)
regio <- do.call(rbind, regio)
t1 <- Sys.time() - t1
saveRDS(regio, "./Figure/2_sim/benchmark/regio.rds")
write.csv(regio, "./Figure/2_sim/benchmark/regio.csv")


### islet
t2 <- Sys.time()
isl <- mclapply(seq(from = 0.1, to = 1, length.out = 30),
                islsimulation,
                MoreArgs = list(
                  n = 1,
                  knot = issim,
                  seurat = surt,
                  senna = seni,
                  giotto = gio,
                  spatial = coord,
                  pval = 0.05),
                mc.preschedule = TRUE,
                mc.cores = 5)
isl <- do.call(rbind, isl)
t2 <- Sys.time() - t2
#saveRDS(isl, "./Figure/benchmark/sim_report/isl.rds")
#write.csv(isl, "./Figure/benchmark/sim_report/isl.csv")


# ver Windows
t2 <- Sys.time()
isl <- lapply(seq(from = 0.1, to = 1, length.out = 30),
              islsimulation,
              n = 1,
              poipar = baseline,
              knot = issim,
              seurat = surt,
              senna = seni,
              giotto = gio,
              sparkarg = spkcoord,
              pval = 0.05)
isl <- do.call(rbind, isl)
t2 <- Sys.time() - t2
saveRDS(isl, "./Figure/2_sim/benchmark/isl.rds")
write.csv(isl, "./Figure/2_sim/benchmark/isl.csv")


# Fig 2. d-e

prog <- readRDS(paste0(path, "dataset/sim/prog.rds"))
regio <- readRDS(paste0(path, "dataset/sim/regio.rds"))
isl <- readRDS(paste0(path, "dataset/sim/isl.rds"))

colnames(prog) <- c("ID", "SENNA", "Seurat", "SPARK-X", "Giotto (k-means)", "Giotto (rank)", "SS")
colnames(regio) <- c("ID", "SENNA", "Seurat", "SPARK-X", "Giotto (k-means)", "Giotto (rank)", "SS")
colnames(isl) <- c("ID", "SENNA", "Seurat", "SPARK-X", "Giotto (k-means)", "Giotto (rank)", "SS")



drsim <- dplyr::filter(prog, ID == "DR_sim") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "DR") %>%
  mutate(ID = "Simulated")

drprm <- dplyr::filter(prog, ID == "DR_prm") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "DR")%>%
  mutate(ID = "Permuted")

drprog <- rbind(drsim, drprm) %>%
  mutate(Study = "Progression") %>%
  filter(between(SS, 0.1, 0.4))
  


drsim <- dplyr::filter(regio, ID == "DR_sim") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "DR") %>%
  mutate(ID = "Simulated")

drprm <- dplyr::filter(regio, ID == "DR_prm") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "DR") %>%
  mutate(ID = "Permuted")

drregio <- rbind(drsim, drprm) %>%
  mutate(Study = "Regionation") %>%
  filter(between(SS, 0.1, 0.4))


drsim <- dplyr::filter(isl, ID == "DR_sim") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "DR") %>%
  mutate(ID = "Simulated")

drprm <- dplyr::filter(isl, ID == "DR_prm") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "DR")%>%
  mutate(ID = "Permuted")

drisl <- rbind(drsim, drprm) %>%
  mutate(Study = "Islet") %>%
  filter(between(SS, 0.1, 0.9))

dr <- rbind(drprog, drregio, drisl)

dr <- dr %>%
  mutate(alg = factor(alg, levels = c("Seurat", 
                                      "SENNA", 
                                      "SPARK-X", 
                                      "Giotto (k-means)", 
                                      "Giotto (rank)")),
         Study = factor(Study, levels = c("Islet",
                                          "Regionation",
                                          "Progression")))


rm("drprm", "drsim", "drprog", "drregio", 'drisl');gc()

pd <- ggplot() +
  geom_line(aes(SS, DR, colour = alg, lty = ID), 
            data = dr, alpha = 0.8) +
  scale_linetype_manual(values = c("Simulated" = "solid", 
                                   "Permuted" = "dashed"),
                        breaks = c("Simulated", "Permuted")) +
  scale_color_manual(values = c("Giotto (k-means)" = "#67b665", 
                                "Seurat" =  "#1b74bc", 
                                "SENNA" = "#ec1b24",
                                "Giotto (rank)" = "#6a3d9a", 
                                "SPARK-X"= "#fa8307"),
                     breaks = c("SENNA",
                                "Seurat",
                                "SPARK-X",
                                "Giotto (k-means)",
                                "Giotto (rank)")) +
  facet_wrap(~Study, ncol = 3, scales = "free_x") +
  labs(x = "Signal Strength",
       y = "Detection Rate",
       lty = "ID",
       color = "Method") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "#dfdfdf"),
        legend.position = "top",
        legend.title = element_text(face = "bold"))

get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

CommonLegend<-get_only_legend(pd)

pd <- ggplot() +
  geom_line(aes(SS, DR, colour = alg, lty = ID), 
            data = dr, alpha = 0.8) +
  scale_linetype_manual(values = c("Simulated" = "solid", 
                                   "Permuted" = "dashed"),
                        breaks = c("Simulated", "Permuted")) +
  scale_color_manual(values = c("Giotto (k-means)" = "#67b665", 
                                "Seurat" =  "#1b74bc", 
                                "SENNA" = "#ec1b24",
                                "Giotto (rank)" = "#6a3d9a", 
                                "SPARK-X"= "#fa8307"),
                     breaks = c("SENNA",
                                "Seurat",
                                "SPARK-X",
                                "Giotto (k-means)",
                                "Giotto (rank)")) +
  facet_wrap(~Study, ncol = 3, scales = "free_x") +
  labs(x = "Signal Strength",
       y = "Detection Rate",
       lty = "ID",
       color = "Method") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "#dfdfdf"),
        legend.position = "none",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(size = 9))



fprp <- dplyr::filter(prog, ID == "FPR_sim") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "FPR") %>%
  mutate(Study = "Progression") %>%
  filter(between(SS, 0.1, 0.4))

fprr <- dplyr::filter(regio, ID == "FPR_sim") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "FPR")%>%
  mutate(Study = "Regionation") %>%
  filter(between(SS, 0.1, 0.4))

fpri <- dplyr::filter(isl, ID == "FPR_sim") %>%
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "alg",
                      values_to = "FPR") %>%
  mutate(Study = "Islet") %>%
  filter(between(SS, 0.1, 0.9))


fprt <- rbind(fprp, fprr, fpri)
fprt <- fprt %>%
  mutate(alg = factor(alg, levels = c("Seurat", 
                                      "SENNA", 
                                      "SPARK-X", 
                                      "Giotto (k-means)", 
                                      "Giotto (rank)")),
         Study = factor(Study, levels = c("Islet",
                                          "Regionation",
                                          "Progression")))

rm("fprp", "fprr", "fpri");gc()

pe <- ggplot() +
  geom_line(aes(SS, FPR, colour = alg), 
            data = fprt,
            lty = "solid",
            alpha = 0.8) +
  scale_color_manual(values = c("Giotto (k-means)" = "#67b665", 
                                "Seurat" =  "#1b74bc", 
                                "SENNA" = "#ec1b24",
                                "Giotto (rank)" = "#6a3d9a", 
                                "SPARK-X"= "#fa8307"),
                     breaks = c("SENNA",
                                "Seurat",
                                "SPARK-X",
                                "Giotto (k-means)",
                                "Giotto (rank)")) +
  facet_wrap(~Study, ncol = 3, scales = "free_x") +
  labs(x = "Signal Strength",
       y = "FPR (1-Specificity)") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "#dfdfdf"),
        legend.position = "none",
        axis.title = element_text(size = 9))


pde <- grid.arrange(pd, pe, CommonLegend, nrow = 3,
                    layout_matrix = rbind(1,1,1,1,1,1,
                                          2,2,2,2,2,2,
                                          3))
ggsave(
  paste0(path, "SENNA/Fig/2_sim/pde_fpr_zoom.tif"),
  pde, 
  dpi = 600,
  width = 18 / 2, 
  height = 10 / 2)


#### SS plot (Results)----

s1 <- Load10X_Spatial(
  data.dir = paste0(path, 
                    "dataset/SNUH/Thymus9/outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Visium",
  filter.matrix = TRUE,
  image = NULL)

if(min(s1$nCount_Spatial) == 0) s1 <- subset(s1, nCount_Spatial > 0)
surt <- s1
s1 <- SCTransform(s1, assay = "Spatial", variable.features.n = 2000, verbose = FALSE)

sen <- SENNA_Visium(s1,
                    slice_name = "Visium")
rm(s1); gc()

AppDat(sen, 
       image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"),
       image_resolution = "lowres")
#knot_picker()



prsim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/svgp.csv"),
                  header = TRUE)
rgsim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/svgr.csv"),
                  header = TRUE)
issim <- read.csv(paste0(path, "dataset/knots/thy/s1/sim/svgi.csv"),
                  header = TRUE)

seni <- TrimmedCurve(sen, issim, type = "islet")
senr <- FullCurve(sen, rgsim, type = "spline")
senp <- FullCurve(sen, prsim, type = "spline")

seni <- GetCurveParam(seni)
senr <- GetCurveParam(senr)
senp <- GetCurveParam(senp)

seni <- TissueRegionation(seni)
senr <- TissueRegionation(senr)

set.seed(1)
baseline <- rnorm(n = 2000, mean = 10, sd = 1)

si <- 1

cprog <- senp
ref <- cprog@Coord[["Spatial"]]
ref[["t"]] <- ref[["t"]] - min(ref[["t"]])
dist <- max(ref[["distance"]]) - ref[["distance"]]
dist <- (dist / max(dist))^2
t <- ref[["t"]] / max(ref[["t"]])

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(t * dist, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
prog_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cprog@Gene[["Spatial"]][1:2000]))))
gs <- rownames(prog_count)[1:1000]
gc <- rownames(prog_count)[1001:2000]

surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(prog_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(prog_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senp <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senp <- FullCurve(senp, prsim, "spline")
senp <- GetCurveParam(senp)

ppd <- senp@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senp@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Progression",
         ID = "SS = 1")
prtj <- SENNA::crvtrjry(senp) %>%
  select(X1, X2) %>%
  mutate(sce = "Progression",
         ID = "SS = 1")



cregio <- senr
ref <- cregio@Coord[["Spatial"]]
csd <- (ref[["distance"]] - min(ref[["distance"]])) / 
  (max(ref[["distance"]]) - min(ref[["distance"]]))
csd[csd == 0] <- 1e-10


grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
regio_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cregio@Gene[["Spatial"]][1:2000]))))

gs <- rownames(regio_count)[1:1000]
gc <- rownames(regio_count)[1001:2000]
surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(regio_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(regio_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senr <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senr <- FullCurve(senr, rgsim, "spline")
senr <- GetCurveParam(senr)

rpd <- senr@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senr@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Regionation",
         ID = "SS = 1")
rgtj <- SENNA::crvtrjry(senr) %>%
  select(X1, X2) %>%
  mutate(sce = "Regionation",
         ID = "SS = 1")


cisl <- seni
ref <- cisl@Coord[["Spatial"]]
csd <- ref[["distance"]]
csd[csd >= 0] <- 0
csd <- -csd
csd <- (csd - min(csd)) / (max(csd) - min(csd))

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref)) 

set.seed(seed = 1)
isl_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cisl@Gene[["Spatial"]][1:2000]))))

gs <- rownames(isl_count)[1:1000]
gc <- rownames(isl_count)[1001:2000]

surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(isl_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(isl_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
seni <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
seni <- TrimmedCurve(seni, issim, "islet")
seni <- GetCurveParam(seni)
seni <- TissueRegionation(seni)

ipd <- seni@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = seni@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Islet",
         ID = "SS = 1")
istj <- SENNA::crvtrjry(seni) %>%
  select(X1, X2) %>%
  mutate(sce = "Islet",
         ID = "SS = 1")

pad <- rbind(ppd,rpd, ipd)
pad <- pad %>%
  mutate(sce = factor(sce, level = c("Islet", "Regionation", "Progression")),
         ID = as.factor(ID))
trj <- rbind(prtj, rgtj, istj)
trj <- trj %>%
  mutate(sce = factor(sce, level = c("Islet", "Regionation", "Progression")),
         ID = as.factor(ID))

pa <- ggplot() +
  geom_point(aes(X1, X2, color = Counts),
             data = pad,
             size = .5,
             alpha = 1) +
  scale_color_viridis_c(
    values = scales::rescale(
      stats::quantile(pad$Counts,
                      probs = c(0, 1)))) +
  geom_point(aes(X1, X2),
             data = trj,
             color = "#000000",
             size = 0.1) +
  theme_light() +
  facet_grid(sce ~ ., switch = "y") +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = c(1, 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "#dfdfdf")) +
  guides(color = "none")
pa
ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/pat_res.tif"),
       pa, 
       dpi = 600,
       width = 7 * .4, 
       height = 9 * .4)





#### SS plot (Supple)----

cprog <- senp
ref <- cprog@Coord[["Spatial"]]
ref[["t"]] <- ref[["t"]] - min(ref[["t"]])
dist <- max(ref[["distance"]]) - ref[["distance"]]
dist <- (dist / max(dist))^2
t <- ref[["t"]] / max(ref[["t"]])

si <- 0.1

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(t * dist, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
prog_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cprog@Gene[["Spatial"]][1:2000]))))
gs <- rownames(prog_count)[1:1000]
gc <- rownames(prog_count)[1001:2000]
surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(prog_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(prog_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senp <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senp <- FullCurve(senp, prsim, "spline")
senp <- GetCurveParam(senp)

ppd <- senp@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senp@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Progression",
         ID = "SS = 0.1")
prtj <- SENNA::crvtrjry(senp) %>%
  select(X1, X2) %>%
  mutate(sce = "Progression",
         ID = "SS = 0.1")

si <- 0.6

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(t * dist, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
prog_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cprog@Gene[["Spatial"]][1:2000]))))

gs <- rownames(prog_count)[1:1000]
gc <- rownames(prog_count)[1001:2000]
surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(prog_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(prog_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senp <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senp <- FullCurve(senp, prsim, "spline")
senp <- GetCurveParam(senp)

ppd0 <- senp@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senp@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Progression",
         ID = "SS = 0.6")
prtj0 <- SENNA::crvtrjry(senp) %>%
  select(X1, X2) %>%
  mutate(sce = "Progression",
         ID = "SS = 0.6")


si <- 1

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(t * dist, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
prog_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    c(gs, gc))))

surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(prog_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(prog_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senp <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senp <- FullCurve(senp, prsim, "spline")
senp <- GetCurveParam(senp)

ppd1 <- senp@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senp@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Progression",
         ID = "SS = 1.0")
prtj1 <- SENNA::crvtrjry(senp) %>%
  select(X1, X2) %>%
  mutate(sce = "Progression",
         ID = "SS = 1.0")

rm("surt_", "prog_count", "grad_index"); gc()


cregio <- sen
cregio <- FullCurve(senna = cregio,
                    knot_df = rgsim,
                    type = "spline")
cregio <- GetCurveParam(cregio)
cregio <- TissueRegionation(cregio)
ref <- cregio@Coord[["Spatial"]]
csd <- (ref[["distance"]] - min(ref[["distance"]])) / 
  (max(ref[["distance"]]) - min(ref[["distance"]]))
csd[csd == 0] <- 1e-10


si <- 0.1

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
regio_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cregio@Gene[["Spatial"]][1:2000]))))

gs <- rownames(regio_count)[1:1000]
gc <- rownames(regio_count)[1001:2000]
surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(regio_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(regio_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senr <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senr <- FullCurve(senr, rgsim, "spline")
senr <- GetCurveParam(senr)

rpd <- senr@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senr@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Regionation",
         ID = "SS = 0.1")
rgtj <- SENNA::crvtrjry(senr) %>%
  select(X1, X2) %>%
  mutate(sce = "Regionation",
         ID = "SS = 0.1")


si <- 0.6

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
regio_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cregio@Gene[["Spatial"]][1:2000]))))

gs <- rownames(regio_count)[1:1000]
gc <- rownames(regio_count)[1001:2000]
surt_ <- surt
surt_ <- subset(surt_, features = c(gs, gc))
surt_@assays[["Spatial"]]@layers$counts <-
  as.sparse(regio_count)
surt_@meta.data$nCount_Spatial <- colSums(regio_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senr <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senr <- FullCurve(senr, rgsim, "spline")
senr <- GetCurveParam(senr)

rpd0 <- senr@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senr@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Regionation",
         ID = "SS = 0.6")
rgtj0 <- SENNA::crvtrjry(senr) %>%
  select(X1, X2) %>%
  mutate(sce = "Regionation",
         ID = "SS = 0.6")


si <- 1

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
regio_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    c(gs, gc))))
surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(regio_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(regio_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
senr <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
senr <- FullCurve(senr, rgsim, "spline")
senr <- GetCurveParam(senr)


rpd1 <- senr@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = senr@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Regionation",
         ID = "SS = 1.0")
rgtj1 <- SENNA::crvtrjry(senr) %>%
  select(X1, X2) %>%
  mutate(sce = "Regionation",
         ID = "SS = 1.0")


rm("surt_", "ref", "csd", "regio_count", "grad_index"); gc()



cisl <- sen
cisl <- TrimmedCurve(senna = cisl,
                     knot_df = issim,
                     type = "islet")
cisl <- GetCurveParam(cisl)
cisl <- TissueRegionation(cisl)
ref <- cisl@Coord[["Spatial"]]
csd <- ref[["distance"]]
csd[csd >= 0] <- 0
csd <- -csd
csd <- (csd - min(csd)) / (max(csd) - min(csd))


si <- 0.1

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref)) 

set.seed(seed = 1)
isl_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cisl@Gene[["Spatial"]][1:2000]))))

gs <- rownames(isl_count)[1:1000]
gc <- rownames(isl_count)[1001:2000]

surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(isl_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(isl_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
seni <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
seni <- TrimmedCurve(seni, issim, "islet")
seni <- GetCurveParam(seni)
seni <- TissueRegionation(seni)

ipd <- seni@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = seni@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Islet",
         ID = "SS = 0.1")
istj <- SENNA::crvtrjry(seni) %>%
  select(X1, X2) %>%
  mutate(sce = "Islet",
         ID = "SS = 0.1")


si <- 0.6

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref)) 

set.seed(seed = 1)
isl_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    colnames(cisl@Gene[["Spatial"]][1:2000]))))

gs <- rownames(isl_count)[1:1000]
gc <- rownames(isl_count)[1001:2000]

surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(isl_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(isl_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
seni <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
seni <- TrimmedCurve(seni, issim, "islet")
seni <- GetCurveParam(seni)
seni <- TissueRegionation(seni)

ipd0 <- seni@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = seni@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Islet",
         ID = "SS = 0.6")
istj0 <- SENNA::crvtrjry(seni) %>%
  select(X1, X2) %>%
  mutate(sce = "Islet",
         ID = "SS = 0.6")



si <- 1

grad_index <- c(si * baseline[1:500],
                - si * baseline[501:1000],
                rep(0, 1e3))
params <- outer(csd, grad_index) + 
  matrix(rep(baseline, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

set.seed(seed = 1)
isl_count <- t(
  matrix(
    rpois(length(params), lambda = params), 
    nrow = nrow(params),
    dimnames = list(rownames(ref),
                    c(gs, gc))))
surt_ <- surt
surt_ <- CreateSeuratObject(
  counts = as.sparse(isl_count),
  assay = "Spatial",
  meta.data = surt@meta.data)
surt_@images <- surt@images
surt_@meta.data$nCount_Spatial <- colSums(isl_count)
surt_ <- NormalizeData(surt_, verbose = FALSE)
surt_ <- ScaleData(surt_, verbose = FALSE)
seni <- SENNA_Visium(surt_,
                     assay = "Spatial",
                     all_genes = TRUE,
                     slice_name = "Visium")
seni <- TrimmedCurve(seni, issim, "islet")
seni <- GetCurveParam(seni)
seni <- TissueRegionation(seni)

ipd1 <- sen@Coord[["Spatial"]] %>%
  select(X1, X2) %>%
  mutate(Counts = seni@Gene[["Spatial"]][[gs[which.max(baseline)]]],
         sce = "Islet",
         ID = "SS = 1.0")
istj1 <- SENNA::crvtrjry(seni) %>%
  select(X1, X2) %>%
  mutate(sce = "Islet",
         ID = "SS = 1.0")



rm("surt_", "ref", "csd", "isl_count", "grad_index"); gc()



pad <- rbind(ppd, ppd0, ppd1, 
             rpd, rpd0, rpd1, 
             ipd, ipd0, ipd1)
pad <- pad %>%
  mutate(sce = factor(sce, level = c("Islet", "Regionation", "Progression")),
         ID = factor(ID, level = c("SS = 0.1", "SS = 0.6", "SS = 1.0")))
trj <- rbind(prtj, prtj0, prtj1, 
             rgtj, rgtj0, rgtj1, 
             istj, istj0, istj1)
trj <- trj %>%
  mutate(sce = factor(sce, level = c("Islet", "Regionation", "Progression")),
         ID = factor(ID, level = c("SS = 0.1", "SS = 0.6", "SS = 1.0")))

pa <- ggplot() +
  geom_point(aes(X1, X2, color = Counts),
             data = pad,
             size = .5,
             alpha = 1) +
  scale_color_viridis_c(
    values = scales::rescale(
      stats::quantile(pad$Counts,
                      probs = c(0, 1)))) +
  labs(color = "Counts (normalized)") +
  guides(color = guide_colorbar(
    label = FALSE,
    barwidth = 5,
    barheight = .7)) +
  geom_point(aes(X1, X2),
             data = trj,
             color = "#000000",
             size = 0.1) +
  theme_light() +
  facet_grid(ID ~ sce, switch = "y") +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = c(1, 1),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "#dfdfdf"))
pa

ggsave(
  paste0(path, "SENNA/Supplementary/sfig/S6_patterns.tif"),
  pa, 
  dpi = 600,
  width = 21, 
  height = 19, 
  units = "cm")








