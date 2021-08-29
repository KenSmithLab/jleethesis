##Nasal - Macrophage Analysis

library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)

setwd("~/Jacinta/Nasal/Pre_processing/")
load("myeloid_semifinal - Nasal.RData")

setwd("~/Jacinta/gene_sets/")
genelist <- list.files(full.names = TRUE, recursive = T, patt ="*.txt")
genelistfinal <- lapply(genelist, read.delim)
setwd("~/Jacinta/Nasal/Macrophage Analysis/")


bio_TGFb <- genelistfinal[[1]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(bio_TGFb), name = "bio_TGFb")
CD16neg_vs_pos_UP <- genelistfinal[[2]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(CD16neg_vs_pos_UP), name = "CD16neg_vs_pos_UP")
M1_vs_M2_DN <- genelistfinal[[3]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(M1_vs_M2_DN), name = "M1_vs_M2_DN")
M1_vs_M2_UP <- genelistfinal[[4]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(M1_vs_M2_UP), name = "M1_vs_M2_UP")
Mono_vs_M1_UP <- genelistfinal[[5]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(Mono_vs_M1_UP), name = "Mono_vs_M1_UP")
Mono_vs_M2_UP <- genelistfinal[[6]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(Mono_vs_M2_UP), name = "Mono_vs_M2_UP")
hall_TGFb <- genelistfinal[[7]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(hall_TGFb), name = "hall_TGFb")

p1 <- FeaturePlot(myeloid.semifinal, features = "bio_TGFb1", order = TRUE, split.by = "timepoint")
p2 <- FeaturePlot(myeloid.semifinal, features = "CD16neg_vs_pos_UP1", order = TRUE)
p3 <- FeaturePlot(myeloid.semifinal, features = "M1_vs_M2_DN1", order = TRUE)
p4 <- FeaturePlot(myeloid.semifinal, features = "M1_vs_M2_UP1", order = TRUE)
p5 <- FeaturePlot(myeloid.semifinal, features = "Mono_vs_M1_UP1", order = TRUE)
p6 <- FeaturePlot(myeloid.semifinal, features = "Mono_vs_M2_UP1", order = TRUE)
p7 <- FeaturePlot(myeloid.semifinal, features = "hall_TGFb1", order = TRUE)
p8 <- p1|p2|p3|p4|p5|p6|p7
p8

#Differential gene expression between Relapse and Remission
Idents(myeloid.semifinal) <- "timepoint"
relapse_markers <- FindMarkers(myeloid.semifinal, ident.1 = "D1", ident.2 = "W12", verbose = FALSE)

tgf_genes <- as.vector(unlist(read.delim("~/Jacinta/gene_sets/TGFb_genes.txt")))
tgf_markers <- relapse_markers[rownames(relapse_markers) %in% tgf_genes,]
write.csv(tgf_markers, "tgfb_markers.csv")

#Differential gene expression between HC and Relapse
hc_markers <- FindMarkers(myeloid.semifinal, ident.1 = "control", ident.2 = "D1", verbose = FALSE)

hc_tgf_markers <- hc_markers[rownames(hc_markers) %in% tgf_genes,]
