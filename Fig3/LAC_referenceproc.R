library(SENNA)
library(Seurat)
library(ggplot2)
library(dplyr)

# Reference pre-processing


load("C:/Users/lmseo/RWD/SPATCH/data/sc/hca_lung_atlas/droplet_normal_lung_seurat_ntiss10x/droplet_normal_lung_blood_seurat_ntiss10x.P2.anno.20191002.RC4.Robj")
ntiss10x.P2.anno <- UpdateSeuratObject(ntiss10x.P2.anno)
common_genes <- intersect(VariableFeatures(object),
                          Features(ntiss10x.P2.anno))
ntiss10x.P2.anno <- subset(ntiss10x.P2.anno, features = common_genes)
saveRDS(ntiss10x.P2.anno, "E:/dataset/Reference/droplet_normal_lung_blood_filtered/vhd1/luad_ref_P2_d1.rds")
rm(ntiss10x.P2.anno); gc()