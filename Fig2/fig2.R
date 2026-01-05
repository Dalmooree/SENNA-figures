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






