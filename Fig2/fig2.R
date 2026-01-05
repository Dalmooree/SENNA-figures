suppressPackageStartupMessages(
  suppressWarnings( {
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(ggsci)
    library(tidyr)
    library(gridExtra)
    library(ggbreak)
    library(forcats)
  }))

path <- getwd()

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
          bg_dot_size = .4,
          bg_dot_alpha = .9,
          line_size = .7,
          line_color = ggsci::pal_npg("nrc")(1),
          knots_color = ggsci::pal_npg("nrc")(1),
          knots_size = 0.1) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4,
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[1, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4, 
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[2, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4, 
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[3, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4, 
            alpha = 0.7,
            data = rbind(qry[,c("X1","X2")], bks[4, c("X1", "X2")])) +
  geom_point(aes(X1, X2), data = bks,
             size = .7, col = ggsci::pal_npg("nrc")(1)) +
  geom_point(aes(X1, X2), data = qry, 
             shape = 21, col = ggsci::pal_npg("nrc")(3)[3], size = .7,
             fill = ggsci::pal_npg("nrc")(3)[3]) +
  geom_point(aes(X1, X2), 
             data = data.frame(X1 = plug_coef(qry$t, sen@CurveAxis$fun$x.coef), X2 = plug_coef(qry$t, sen@CurveAxis$fun$y.coef)),
             shape = 21, col = ggsci::pal_npg("nrc")(3)[3], size = .7,
             fill = ggsci::pal_npg("nrc")(3)[3]) +
  geom_line(aes(X1, X2), linetype = "solid", 
            col = ggsci::pal_npg("nrc")(3)[3], 
            linewidth = .5, 
            data = rbind(qry[,c("X1","X2")], data.frame(X1 = plug_coef(qry$t, sen@CurveAxis$fun$x.coef), X2 = plug_coef(qry$t, sen@CurveAxis$fun$y.coef)))) +
  ggrepel::geom_label_repel(aes(X1, X2, 
                               label = paste0("t=", round(t, 1))), 
                           data = bks,
                           box.padding = 0.1,
                           size = 2.5,
                           fill = alpha("#ffffff", 0.5),
                           color = "#000000",
                           seed = 123) + 
  ggrepel::geom_label_repel(aes(X1, X2, 
                               label = paste0("t=", round(t, 1))), 
                           data = data.frame(X1 = plug_coef(qry$t, sen@CurveAxis$fun$x.coef), 
                                             X2 = plug_coef(qry$t, sen@CurveAxis$fun$y.coef),
                                             t = sen@Coord$Spatial$t[rownames(sen@Coord$Spatial) %in% rownames(qry)]),
                           fontface = "bold",
                           box.padding = 0.1,
                           nudge_y = 0.02,
                           fill = alpha("#ffffff", 0.8),
                           color = "#000000",,
                           label.padding = unit(0.3, "lines"),
                           size = 3,
                           seed = 1) +
  labs(x = "X", y = "Y") + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/prog_ca2.tif"),
       dpi = 600, width = 11 / 4, height = 0.7 * 11 / 4)


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
             col = "#444444", alpha = 0.7,
             linewidth = .2, linetype = "longdash") +
  geom_vline(xintercept = bks$t[2],
             col = "#444444", alpha = 0.7,
             linewidth = .2, linetype = "longdash") +
  geom_vline(xintercept = bks$t[3],
             col = "#444444", alpha = 0.7,
             linewidth = .2, linetype = "longdash") +
  geom_vline(xintercept = bks$t[4],
             col = "#444444", alpha = 0.7,
             linewidth = .2, linetype = "longdash") +
  geom_vline(xintercept = qry$t,
             col = ggsci::pal_npg("nrc")(3)[3],
             linewidth = .5, linetype = "longdash") +
  geom_line(aes(t, Distance), data = td,
            col = ggsci::pal_npg("nrc")(4)[4], linewidth = .7) +
  labs(x = "Curve Parameter (t)") +
  scale_x_continuous(breaks = c(0, bks$t[4]),
                     labels = c("0",
                                format(bks$t[4], digits = 2))) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  ggrepel::geom_label_repel(
    aes(x = t, y = distance, label = paste0("t=", round(t, 1))),
    data = bks,
    box.padding = 0.1,
    size = 2.5,
    fill = alpha("#ffffff", 0.8),
    color = "#000000",
    direction = "y",
    seed = 111
  ) +
  geom_label(
    aes(x = t, y = .5, label = paste0("t=", round(t, 1))),
    data = qry,
    fontface = "bold",
    fill = alpha("#ffffff", 0.8),
    color = "#000000",,
    label.padding = unit(0.3, "lines"),
    size = 3
  )

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/prog_dist2.tif"),
       dpi = 600, width = 11/4, height = 0.7*11/4)


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
          bg_dot_size = .4,
          bg_dot_alpha = .9,
          line_size = .7,
          line_color = ggsci::pal_npg("nrc")(1),
          knots_color = ggsci::pal_npg("nrc")(1),
          knots_size = 0.1) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[1, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[2, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[3, c("X1", "X2")])) +
  geom_line(aes(X1, X2), linetype = "longdash", 
            col = "#444444", 
            linewidth = .4, alpha = .7, 
            data = rbind(qry[,c("X1","X2")], bks[4, c("X1", "X2")])) +
  geom_point(aes(X1, X2), data = bks,
             size = .7, col = ggsci::pal_npg("nrc")(1)) +
  geom_point(aes(X1, X2), data = qry, 
             shape = 21, fill = ggsci::pal_npg("nrc")(3)[3], 
             col = ggsci::pal_npg("nrc")(3)[3], size = 1.7) +
  geom_point(aes(X1, X2), 
             data = data.frame(X1 = trplug_coef(qry$t, senc@CurveAxis$fun$x.coef), X2 = trplug_coef(qry$t, senc@CurveAxis$fun$y.coef)),
             shape = 21, fill = ggsci::pal_npg("nrc")(3)[3], size = .7,
             col = ggsci::pal_npg("nrc")(3)[3]) +
  geom_line(aes(X1, X2), linetype = "solid", 
            col = ggsci::pal_npg("nrc")(3)[3], 
            linewidth = .5, 
            data = rbind(qry[,c("X1","X2")], 
                         data.frame(X1 = trplug_coef(qry$t, senc@CurveAxis$fun$x.coef), 
                                    X2 = trplug_coef(qry$t, senc@CurveAxis$fun$y.coef)))) +
  ggrepel::geom_label_repel(aes(X1, X2, 
                               label = paste0("t=", round(t, 1))), 
                           data = bks,,
                           point.padding = 0.01,
                           box.padding = 0.1,
                           size = 2.5,
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
                           label.padding = unit(0.3, "lines"),
                           size = 3,
                           seed = 1) +
  labs(x = "X", y = "Y")+ 
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/isl_ca1.tif"),
       dpi = 600, width = 11/4, height = 0.7*11/4)

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
             col = "#444444", alpha = 0.7,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[2],
             col = "#444444", alpha = 0.7,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[3],
             col = "#444444", alpha = 0.7,
             linewidth = .3, linetype = "longdash") +
  geom_vline(xintercept = bks$t[4],
             col = "#444444", alpha = 0.7,
             linewidth = .3, linetype = "longdash") +
  geom_line(aes(t, Distance), data = td,
            col = ggsci::pal_npg("nrc")(4)[4], linewidth = .7) +
  geom_vline(xintercept = qry$t,
             col = ggsci::pal_npg("nrc")(3)[3],
             linewidth = .5, linetype = "longdash")+ 
  labs(x = "Curve Parameter (t)") +
  scale_x_continuous(breaks = c(1, max(td)),
                     labels = c("1", 
                                format(max(td), digits = 1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  ggrepel::geom_label_repel(
    aes(x = t, y = distance, label = paste0("t=", round(t, 1))),
    data = bks,
    box.padding = 0.1,
    size = 2.5,
    fill = alpha("#ffffff", 0.8),
    color = "#000000",
    direction = "y",
    seed = 123
  ) +
  geom_label(
    aes(x = t, y = .23, label = paste0("t=", round(t, 2))),
    data = qry,
    fontface = "bold",
    fill = alpha("#ffffff", 0.8),
    color = "#000000",,
    label.padding = unit(0.3, "lines"),
    size = 3
  )

ggsave(paste0(path,
              "SENNA/Fig/2_sim/simulation_benchmarking/expr1/isl_dist2.tif"),
       dpi = 600, width = 11/4, height = 0.7*11/4)


rm(list = ls());gc()
## Fig 2, d-e

# Load data ------------------------------------------------
si <- seq(from = 0.1, to = 1, length.out = 30)
dpath <- "data/bhmk.res"

# Initialize empty data frames
pt <- data.frame()
p.time <- data.frame()

## Prog - Seurat
for(siv in si){
  siv <- round(siv, 2)
  tmp <- readRDS(file.path(dpath, "p", "SEURAT", paste0("res_seurat_SI", siv, ".rds")))
  
  pt <- rbind(pt, tmp %>%
                select(ID, SI, value = SEURAT) %>%
                mutate(method = "SEURAT"))
  
  p.time <- rbind(p.time, tmp %>%
                    filter(ID == "Power_sim") %>%
                    select(SI, total_time, mapping_time) %>%
                    mutate(method = "Seurat"))
}

# Prog - Other R methods
rmtds <- c("SPARK", "SPARK_X", "nnSVG", "SENNA")
for(m in rmtds){
  f <- if(m == "SPARK_X") "sparkx" else tolower(m)
  d <- if(m == "SPARK_X") "SPARKX" else toupper(m)
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- readRDS(file.path(dpath, "p", d, paste0("res_", f, "_SI", siv, ".rds")))
    
    pt <- rbind(pt, tmp %>%
                  select(ID, SI, value = all_of(toupper(m))) %>%
                  mutate(method = toupper(m)))
    
    p.time <- rbind(p.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = m))
  }
}

# Prog - Giotto (kmeans와 rank 두 가지)
giotto_types <- c("giottok", "giottor")
giotto_cols <- c("GIOTTO_kmeans", "GIOTTO_rank")

for(j in 1:length(giotto_types)){
  gtype <- giotto_types[j]
  gcol <- giotto_cols[j]
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- readRDS(file.path(dpath, "p", "GIOTTO", paste0("res_", gtype, "_SI", siv, ".rds")))
    
    pt <- rbind(pt, tmp %>%
                  select(ID, SI, value = all_of(gcol)) %>%
                  mutate(method = gcol))
    
    p.time <- rbind(p.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = paste0("Giotto_", ifelse(gtype == "giottok", "kmeans", "rank"))))
  }
}

# Prog - Python methods
pymtds <- c("SOMDE", "SPATIALDE2")

for(m in pymtds){
  # SpatialDE2는 파일명이 spatialde
  f <- if(m == "SPATIALDE2") "spatialde" else tolower(m)
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- read.csv(file.path(dpath, "p", m, paste0("res_", f, "_SI", siv, ".csv")))
    
    pt <- rbind(pt, tmp %>%
                  select(ID, SI, value = all_of(m)) %>%
                  mutate(method = m))
    
    p.time <- rbind(p.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = m))
  }
}

## Regio (r)
rt <- data.frame()
r.time <- data.frame()

# Regio - Seurat
for(siv in si){
  siv <- round(siv, 2)
  tmp <- readRDS(file.path(dpath, "r", "SEURAT", paste0("res_seurat_SI", siv, ".rds")))
  
  rt <- rbind(rt, tmp %>%
                select(ID, SI, value = SEURAT) %>%
                mutate(method = "SEURAT"))
  
  r.time <- rbind(r.time, tmp %>%
                    filter(ID == "Power_sim") %>%
                    select(SI, total_time, mapping_time) %>%
                    mutate(method = "Seurat"))
}

# Regio - Other R methods
for(m in rmtds){
  f <- if(m == "SPARK_X") "sparkx" else tolower(m)
  d <- if(m == "SPARK_X") "SPARKX" else toupper(m)
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- readRDS(file.path(dpath, "r", d, paste0("res_", f, "_SI", siv, ".rds")))
    
    rt <- rbind(rt, tmp %>%
                  select(ID, SI, value = all_of(toupper(m))) %>%
                  mutate(method = toupper(m)))
    
    r.time <- rbind(r.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = m))
  }
}

# Regio - Giotto
for(j in 1:length(giotto_types)){
  gtype <- giotto_types[j]
  gcol <- giotto_cols[j]
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- readRDS(file.path(dpath, "r", "GIOTTO", paste0("res_", gtype, "_SI", siv, ".rds")))
    
    rt <- rbind(rt, tmp %>%
                  select(ID, SI, value = all_of(gcol)) %>%
                  mutate(method = gcol))
    
    r.time <- rbind(r.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = paste0("Giotto_", ifelse(gtype == "giottok", "kmeans", "rank"))))
  }
}

# Regio - Python methods
for(m in pymtds){
  f <- if(m == "SPATIALDE2") "spatialde" else tolower(m)
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- read.csv(file.path(dpath, "r", m, paste0("res_", f, "_SI", siv, ".csv")))
    
    rt <- rbind(rt, tmp %>%
                  select(ID, SI, value = all_of(m)) %>%
                  mutate(method = m))
    
    r.time <- rbind(r.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = m))
  }
}

## Islet (i)
it <- data.frame()
i.time <- data.frame()

# Islet - Seurat
for(siv in si){
  siv <- round(siv, 2)
  tmp <- readRDS(file.path(dpath, "i", "SEURAT", paste0("res_seurat_SI", siv, ".rds")))
  
  it <- rbind(it, tmp %>%
                select(ID, SI, value = SEURAT) %>%
                mutate(method = "SEURAT"))
  
  i.time <- rbind(i.time, tmp %>%
                    filter(ID == "Power_sim") %>%
                    select(SI, total_time, mapping_time) %>%
                    mutate(method = "Seurat"))
}

# Islet - Other R methods
for(m in rmtds){
  f <- if(m == "SPARK_X") "sparkx" else tolower(m)
  d <- if(m == "SPARK_X") "SPARKX" else toupper(m)
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- readRDS(file.path(dpath, "i", d, paste0("res_", f, "_SI", siv, ".rds")))
    
    it <- rbind(it, tmp %>%
                  select(ID, SI, value = all_of(toupper(m))) %>%
                  mutate(method = toupper(m)))
    
    i.time <- rbind(i.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = m))
  }
}

# Islet - Giotto
for(j in 1:length(giotto_types)){
  gtype <- giotto_types[j]
  gcol <- giotto_cols[j]
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- readRDS(file.path(dpath, "i", "GIOTTO", paste0("res_", gtype, "_SI", siv, ".rds")))
    
    it <- rbind(it, tmp %>%
                  select(ID, SI, value = all_of(gcol)) %>%
                  mutate(method = gcol))
    
    i.time <- rbind(i.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = paste0("Giotto_", ifelse(gtype == "giottok", "kmeans", "rank"))))
  }
}

# Islet - Python methods
for(m in pymtds){
  f <- if(m == "SPATIALDE2") "spatialde" else tolower(m)
  
  for(siv in si){
    siv <- round(siv, 2)
    tmp <- read.csv(file.path(dpath, "i", m, paste0("res_", f, "_SI", siv, ".csv")))
    
    it <- rbind(it, tmp %>%
                  select(ID, SI, value = all_of(m)) %>%
                  mutate(method = m))
    
    i.time <- rbind(i.time, tmp %>%
                      filter(ID == "Power_sim") %>%
                      select(SI, total_time, mapping_time) %>%
                      mutate(method = m))
  }
}

# Save data ------------------------------------------------
write.csv(pt, file.path(dpath, "pt.csv"), row.names = FALSE)
write.csv(p.time, file.path(dpath, "p.time.csv"), row.names = FALSE)

write.csv(rt, file.path(dpath, "rt.csv"), row.names = FALSE)
write.csv(r.time, file.path(dpath, "r.time.csv"), row.names = FALSE)

write.csv(it, file.path(dpath, "it.csv"), row.names = FALSE)
write.csv(i.time, file.path(dpath, "i.time.csv"), row.names = FALSE)



# Rstudio; Load Dataset -----
dpath <- file.path("Fig", "2_sim", "benchmark", "sim_report_new")
pt <- read.csv(file.path(dpath, "pt.csv"))
p.time <- read.csv(file.path(dpath, "p.time.csv"))
 
rt <- read.csv(file.path(dpath, "rt.csv"))
r.time <- read.csv(file.path(dpath, "r.time.csv"))
 
it <- read.csv(file.path(dpath, "it.csv"))
i.time <- read.csv(file.path(dpath, "i.time.csv"))


# DR, FPR -----

pt <- mutate(pt, Study = "Progression")
rt <- mutate(rt, Study = "Regionation")
it <- mutate(it, Study = "Islet")

t <- split(pt, pt$ID %in% c("Power_sim", "Power_prm"))
prp <- t$`TRUE`; prf <- t$`FALSE`; rm(pt)
t <- split(rt, rt$ID %in% c("Power_sim", "Power_prm"))
rrp <- t$`TRUE`; rrf <- t$`FALSE`; rm(rt)
t <- split(it, it$ID %in% c("Power_sim", "Power_prm"))
irp <- t$`TRUE`; irf <- t$`FALSE`; rm(it); rm(t)

dr <- rbind(prp, rrp, irp)
fpr <- rbind(prf, rrf, irf)
rm("prp", "rrp", "irp", "prf", "rrf", "irf"); gc()

dr <- mutate(dr, ID = ifelse(dr$ID == "Power_sim", "Simulated", "Permuted"))
fpr <- mutate(fpr, ID = ifelse(fpr$ID == "FPR_sim", "Simulated", "Permuted"))

dr <- dr %>%
  mutate(method = case_when(
    method == "SEURAT" ~ "Seurat",
    method == "SPARK_X" ~ "SPARK-X",
    method == "NNSVG" ~ "nnSVG",
    method == "SPATIALDE2" ~ "SpatialDE2",
    method == "GIOTTO_kmeans" ~ "Giotto (kmeans)",
    method == "GIOTTO_rank" ~ "Giotto (rank)",
    TRUE ~ method
  ))

fpr <- fpr %>%
  mutate(method = case_when(
    method == "SEURAT" ~ "Seurat",
    method == "SPARK_X" ~ "SPARK-X",
    method == "NNSVG" ~ "nnSVG",
    method == "SPATIALDE2" ~ "SpatialDE2",
    method == "GIOTTO_kmeans" ~ "Giotto (kmeans)",
    method == "GIOTTO_rank" ~ "Giotto (rank)",
    TRUE ~ method
  ))


pals <- ggsci::pal_npg("nrc")(10)[c(8, 2:7, 9, 10)]
names(pals) <- c("SENNA", unique(dr$method)[unique(dr$method) != "SENNA"])


dr <- dr %>%
  mutate(method = factor(method, levels = c("SENNA", unique(dr$method)[unique(dr$method) != "SENNA"])),
         Study = factor(Study, levels = c("Islet",
                                          "Regionation",
                                          "Progression")))

fpr <- fpr %>%
  mutate(method = factor(method, levels = c("SENNA", unique(fpr$method)[unique(fpr$method) != "SENNA"])),
         Study = factor(Study, levels = c("Islet",
                                          "Regionation",
                                          "Progression")))


pd <- ggplot() +
  geom_line(aes(SI, value, colour = method, lty = ID), 
            data = dr, alpha = 0.6) +
  scale_linetype_manual(values = c("Simulated" = "solid", 
                                   "Permuted" = "dotted"),
                        breaks = c("Simulated", "Permuted")) +
  scale_color_manual(values = pals) +
  facet_wrap(~Study, ncol = 3, scales = "free_x") +
  labs(x = "Signal Strength",
       y = "Detection Rate",
       lty = "ID",
       color = "Method") +
  theme_light() +
  theme(strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "#dfdfdf"),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 8),
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 7),
        legend.spacing.y = unit(.1, "lines"),
        legend.key.spacing.y = unit(.1, "lines")) +
  guides(lty = guide_legend(ncol = 1))

get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

CommonLegend <- get_only_legend(pd)

pd <- ggplot() +
  geom_line(aes(SI, value, colour = method, lty = ID), 
            data = dr, alpha = 0.6) +
  scale_linetype_manual(values = c("Simulated" = "solid", 
                                   "Permuted" = "dotted"),
                        breaks = c("Simulated", "Permuted")) +
  scale_color_manual(values = pals) +
  facet_wrap(~Study, ncol = 3, scales = "free_x") +
  labs(x = "Signal Strength",
       y = "Detection Rate",
       lty = "ID",
       color = "Method") +
  theme_light() +
  theme(strip.text = element_text(color = "black", size = 7,
                                  margin = margin(t = 2, b = 2)),
        strip.background = element_rect(fill = "#dfdfdf"),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_blank())

pe <- ggplot() +
  geom_line(aes(SI, value, colour = method, lty = ID), 
            data = fpr, alpha = 0.8) +
  scale_linetype_manual(values = c("Simulated" = "solid", 
                                   "Permuted" = "dashed"),
                        breaks = c("Simulated", "Permuted")) +
  scale_color_manual(values = pals) +
  facet_wrap(~Study, ncol = 3, scales = "free_x") +
  labs(x = "Signal Strength",
       y = "False Positive Rate",
       lty = "ID",
       color = "Method") +
  theme_light() +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

pde <- grid.arrange(pd, pe, CommonLegend, nrow = 3,
                    layout_matrix = rbind(1,1,1,1,1,1,
                                          2,2,2,2,2,2,
                                          3,
                                          3)); pde
  
ggsave(
  plot = pde, 
  filename = file.path(dpath, "res", "pde.tif"),
  dpi = 600,
  width = 12 / 2, 
  height = 7 / 2,
  bg = "white")






# Runtime ----

p.time <- mutate(p.time, Study = "Progression")
r.time <- mutate(r.time, Study = "Regionation")
i.time <- mutate(i.time, Study = "Islet")
tt <- rbind(p.time, r.time, i.time)
rm("p.time", "i.time", "r.time"); gc()

  
tt <- tt %>%
  mutate(method = case_when(
    method == "SPARK_X" ~ "SPARK-X",
    method == "SPATIALDE2" ~ "SpatialDE2",
    method == "Giotto_kmeans" ~ "Giotto (kmeans)",
    method == "Giotto_rank" ~ "Giotto (rank)",
    TRUE ~ method
  ))


tt <- tt %>%
  mutate(method = fct_reorder(method, total_time, .fun = mean, .desc = TRUE)) %>%
  mutate(method = fct_relevel(method, "SENNA", after = Inf))


stt <- tt %>%
  group_by(method) %>%
  summarise(mu = mean(total_time, na.rm = TRUE),
            sigma = sd(total_time, na.rm = TRUE)) %>%
  ungroup()

prt <- ggplot(aes(mu, method), data = stt) +
  geom_col(aes(fill = method)) +
  scale_fill_manual(values = pals,
                    breaks = stt[["method"]]) + 
  geom_errorbar(aes(xmin = mu - sigma,
                    xmax = mu + sigma),
                width = .2) +
  geom_point(aes(x = total_time), 
             data = tt, 
             size = .1, 
             alpha = 0.2, 
             position = position_jitter(height = 0.1)) +
  labs(y = "Methods",
       x = "Runtime (second)") +
  scale_x_continuous(breaks = seq(from = 0, to = 1.9e4, by = 1e3),
                     labels = scales::label_scientific(digits = 2)) +
  scale_x_break(c(2000, 17500)) + 
  theme_light() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7, angle = 30),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8)) +
  guides(fill = "none"); prt

ggsave(
  plot = prt, 
  filename = file.path(dpath, "res", "prt.tif"),
  dpi = 600,
  width = 6 / 2,
  height = 7 / 2,
  bg = "white")


tt <- tt %>%
  mutate(method = fct_reorder(method, total_time, .fun = mean, .desc = FALSE)) %>%
  mutate(method = fct_relevel(method, "SENNA", after = 0L))


stt <- tt %>%
  group_by(method) %>%
  summarise(mu = mean(total_time, na.rm = TRUE),
            sigma = sd(total_time, na.rm = TRUE)) %>%
  ungroup()


prt_sup <- ggplot(aes(method, log2(mu)), data = stt) +
  geom_col(aes(fill = method)) +
  scale_fill_manual(values = pals,
                    breaks = stt[["method"]]) + 
  geom_errorbar(aes(ymin = log2(mu - sigma),
                    ymax = log2(mu + sigma)),
                width = .2) +
  geom_point(aes(y = log2(total_time)), 
             data = tt, 
             size = .3, 
             alpha = 0.2, 
             position = position_jitter(width = 0.1)) +
  labs(x = "Methods",
       y = "Runtime (second, log2-scaled)") +
  theme_light() +
  guides(fill = "none"); prt; prt_sup

ggsave(
  plot = prt_sup, 
  filename = file.path(dpath, "res", "prt_sup.tif"),
  dpi = 600,
  width = 18 / 2,
  height = 7 / 2,
  bg = "white")












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







