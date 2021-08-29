## Nasal - Cluster annotations

library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(tidyverse)
library(rlist)

setwd("~/Jacinta/Nasal/Pre_processing/")

#Identifying broad clusters
load("VascSC_int_clustered - Nasal.RData")
DefaultAssay(VascSC.int.clustered) <- "RNA"

all_markers <-FindAllMarkers(VascSC.int.clustered, 
                             min.pct =  0.25, 
                             min.diff.pct = 0.25)
write.csv(all_markers, "All_markers - Nasal.csv")

top10markers <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
write.csv(top10markers, "Top10markers - Nasal.csv")

pdf("Violin_Marker_plot1.pdf", width = 15, height = 13)
VlnPlot(VascSC.int.clustered, features = c("IL7R", "CCR7", "S100A4", "CD14", "LYZ", 
                                           "MS4A1", "CD8A", "FCGR3A", "MS4A7"))
dev.off()

pdf("Violin_Marker_plot2.pdf", width = 15, height = 13)
VlnPlot(VascSC.int.clustered, features = c("GNLY", "NKG7", "FCER1A", "CST3", "PPBP", 
                                           "S100A12", "SELL", "CD3E", "PRF1"))
dev.off()

pdf("Violin_Marker_plot3.pdf", width = 15, height = 13)
VlnPlot(VascSC.int.clustered, features = c("GPR183", "GZMA", "FOXP3", "IL2RA", 
                                           "SLC4A10", "TRAV1-2", "TRGV9", "TRDV2", 
                                           "KLRF1"))
dev.off()

pdf("Nasal - Canonical_markers.pdf", width = 15, height = 13)
FeaturePlot(VascSC.int.clustered, features = c("CD8A", "CD3E", "MS4A1", "GNLY", "NKG7", "EPCAM",
                                               "CD14", "LYZ", "FCGR3A", "FCER1A", "PPBP", "HBB"), 
            order = TRUE)
dev.off()

#SingleR to tease out general population
library(celldex)
library(SingleR)

hpca <- celldex::HumanPrimaryCellAtlasData()

hpcaResults <- SingleR(test = as.SingleCellExperiment(VascSC.int.clustered), ref = hpca, labels = hpca$label.main)
write.csv(hpcaResults, "Nasal - hcpaResults.csv")

results.test <- read.csv('Nasal - hcpaResults.csv', header = T)
rownames(results.test) = results.test$X
VascSC.int.clustered <- AddMetaData(VascSC.int.clustered, results.test)

DimPlot(VascSC.int.clustered, group.by = 'pruned.labels', label = T, label.size = 5, repel = TRUE, raster = FALSE)

table(VascSC.int.clustered$pruned.labels)


#Renaming clusters
new.cluster.ids <- c("T", "B", "B", "T", "T", "T", "Doublet", "T", "Myeloid",
                     "B", "Epithelial", "Epithelial", "B", "NK", "B", "Myeloid", "Myeloid")

names(new.cluster.ids) <- levels(VascSC.int.clustered)
VascSC.int.clustered <- RenameIdents(VascSC.int.clustered, new.cluster.ids)
table(VascSC.int.clustered@active.ident)

pdf("Nasal - Cell subsets.pdf", width = 6, height = 4)
DimPlot(VascSC.int.clustered, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

#No. of cells per cell type
table(VascSC.int.clustered@active.ident)

#Extract Myeloid Cells
myeloid.subset <- subset(x = VascSC.int.clustered, idents = "Myeloid")
table(myeloid.subset$sample)

save(myeloid.subset, file = "Nasal_myeloid_subset.RData")

DefaultAssay(myeloid.subset) <- "RNA"

#Renormalise samples and integrate
myeloidList <- SplitObject(myeloid.subset, split.by = "sample")

trimMyeloidList <- list.remove(myeloidList, range = c(9, 21, 23))

trimMyeloidList <- lapply(X = trimMyeloidList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = trimMyeloidList)
trimMyeloidList <- lapply(X = trimMyeloidList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 8)
})

myeloid.anchors <- FindIntegrationAnchors(object.list = trimMyeloidList,
                                          anchor.features = features)
save(myeloid.anchors, file = "Myeloid_anchors - Nasal.RData")

load("Myeloid_anchors - Nasal.RData")
myeloid.integrated <- IntegrateData(anchorset = myeloid.anchors, k.weight = 30)
save(myeloid.integrated, file = "Myeloid_integrated - Nasal.RData")

DefaultAssay(myeloid.integrated) <- "integrated"

myeloid.integrated <- ScaleData(myeloid.integrated, verbose = FALSE)
myeloid.integrated <- RunPCA(myeloid.integrated, npcs = 30, verbose = FALSE)

#Determine how many PCs to use, select around inflection point
ElbowPlot(myeloid.integrated)

myeloid.integrated <- RunUMAP(myeloid.integrated, reduction = "pca", dims = 1:17)
myeloid.integrated <- FindNeighbors(myeloid.integrated, reduction = "pca", dims = 1:17)
myeloid.clustered <- FindClusters(myeloid.integrated, resolution = 0.3)

save(myeloid.clustered, file = "Myeloid_clustered - Nasal.RData")

#Re-Annotate
DefaultAssay(myeloid.clustered) <- "RNA"

pdf("Nasal - Myeloid Canonical_markers.pdf", width = 15, height = 13)
FeaturePlot(myeloid.clustered, features = c("CD8A", "CD3E", "MS4A1", "GNLY", "NKG7", "KLRF1",
                                            "CD14", "LYZ", "CD1C", "FCER1A", "PPBP", "HBB"), 
            order = TRUE)
dev.off()

##Plot UMAP
pdf(file = "Nasal - Initial Myeloid UMAP.pdf", width = 10, height = 8)
DimPlot(myeloid.clustered, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

pdf(file = "Nasal - Initial Myeloid UMAP - Timepoint.pdf", width = 15, height = 5)
DimPlot(myeloid.clustered, reduction = "umap", group.by = "timepoint", raster = FALSE)
dev.off()

pdf(file = "Nasal - Initial Myeloid UMAP - Batch.pdf", width = 15, height = 5)
DimPlot(myeloid.clustered, reduction = "umap", group.by = "batch", raster = FALSE)
dev.off()


#SingleR to tease out myeloid population
library(celldex)
library(SingleR)

monaco <- celldex::MonacoImmuneData()

monacoResults <- SingleR(test = as.SingleCellExperiment(myeloid.clustered), ref = monaco, labels = monaco$label.fine)
write.csv(monacoResults, "monacoResults.csv")

results.test <- read.csv('monacoResults.csv', header = T)
rownames(results.test) = results.test$X
myeloid.clustered <- AddMetaData(myeloid.clustered, results.test)

DimPlot(myeloid.clustered, group.by = 'pruned.labels', label = T, label.size = 5, repel = TRUE, raster = FALSE)

table(myeloid.clustered$pruned.labels)

#Subsetting to clear doublets, and other cells then re-normalising and integration
myeloid.trim <- subset(x = myeloid.clustered, idents = c(0, 4, 6, 8))
nasal.sampleno <- table(myeloid.trim$sample)
write.csv(nasal.sampleno, "nasal_myeloid_numbers.csv")

DefaultAssay(myeloid.trim) <- "RNA"

#Renormalise samples and integrate, remove samples with less than 30 cells
#c20_nasal_con, c21_nasal_con, cv01_nasal_d1, cv02_nasal_w12_cite
myeloidtrimList <- SplitObject(myeloid.trim, split.by = "sample")
myeloidtrimList2 <- list.remove(myeloidtrimList, range = c(3, 4, 6, 8, 17, 18))

myeloidtrimList2 <- lapply(X = myeloidtrimList2, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = myeloidtrimList2)
myeloidtrimList2 <- lapply(X = myeloidtrimList2, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 29)
})

myeloid.trim.anchors <- FindIntegrationAnchors(object.list = myeloidtrimList2,
                                               anchor.features = features)
save(myeloid.trim.anchors, file = "Myeloid_trim_anchors - Nasal.RData")

load("Myeloid_trim_anchors - Nasal.RData")
myeloid.trim.integrated <- IntegrateData(anchorset = myeloid.trim.anchors, k.weight = 30)
save(myeloid.trim.integrated, file = "Myeloid_trim_integrated - Nasal.RData")

DefaultAssay(myeloid.trim.integrated) <- "integrated"

myeloid.trim.integrated <- ScaleData(myeloid.trim.integrated, verbose = FALSE)
myeloid.trim.integrated <- RunPCA(myeloid.trim.integrated, npcs = 30, verbose = FALSE)

#Determine how many PCs to use, select around inflection point
ElbowPlot(myeloid.trim.integrated)

myeloid.trim.integrated <- RunUMAP(myeloid.trim.integrated, reduction = "pca", dims = 1:19)
myeloid.trim.integrated <- FindNeighbors(myeloid.trim.integrated, reduction = "pca", dims = 1:19)
save(myeloid.trim.integrated, file = "Myeloid_trim_int_before_cluster - Nasal.RData")
myeloid.trim.clustered <- FindClusters(myeloid.trim.integrated, resolution = 0.5)

save(myeloid.trim.clustered, file = "Myeloid_trim_clustered - Nasal.RData")

#Re-Annotate
load("Myeloid_trim_clustered - Nasal.RData")
DefaultAssay(myeloid.trim.clustered) <- "RNA"

pdf("Nasal - Myeloid_trim_Canonical_markers.pdf", width = 15, height = 13)
FeaturePlot(myeloid.trim.clustered, features = c("CD8A", "CD3E", "MS4A1", "GNLY", "NKG7", "KLRF1",
                                                 "CD14", "LYZ", "CD1C", "FCER1A", "PPBP", "HBB"), 
            order = TRUE)
dev.off()

pdf(file = "Nasal - Myeloid_trim.pdf", width = 5, height = 4)
DimPlot(myeloid.trim.clustered, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

#EnrichR
hpcaResults2 <- SingleR(test = as.SingleCellExperiment(myeloid.trim.clustered), ref = hpca, labels = hpca$label.main)
write.csv(hpcaResults2, "Nasal - hcpaResults2.csv")

results.test2 <- read.csv('Nasal - hcpaResults2.csv', header = T)
rownames(results.test2) = results.test2$X
myeloid.trim.clustered <- AddMetaData(myeloid.trim.clustered, results.test2)

DimPlot(myeloid.trim.clustered, group.by = 'pruned.labels', label = T, label.size = 5, repel = TRUE, raster = FALSE)

table(myeloid.trim.clustered$pruned.labels)
table(myeloid.trim.clustered$sample)



myeloid.semifinal <- myeloid.trim.clustered[,grep("Monocyte|Macrophage|DC", myeloid.trim.clustered$pruned.labels)]

pdf("Nasal - UMAP split by SingleR labels - hpca.pdf", width = 15, height = 5)
DimPlot(myeloid.semifinal, split.by = 'pruned.labels', label = T, label.size = 5, repel = TRUE, raster = FALSE)
dev.off()

pdf("Nasal - UMAP split by timepoint.pdf", width = 15, height = 5)
DimPlot(myeloid.semifinal, split.by = 'timepoint', label = T, label.size = 5, repel = TRUE, raster = FALSE)
dev.off()


table(myeloid.semifinal$sample)
table(myeloid.semifinal@active.ident)

save(myeloid.semifinal, file = "myeloid_semifinal - Nasal.RData")
saveRDS(myeloid.semifinal, file = "myeloid_semifinal.rds")

#Save individual RDS
demo <- read.csv("~/Jacinta/single_cell_phenodata.csv")
filename <- demo$filename
filename <- filename[grep("nasal", filename)]
filename <- filename[-c(3, 4, 6, 8, 9, 12, 19, 21, 23)]

myeloid.list <- SplitObject(myeloid.semifinal, split.by = "sample")
for (i in 1:length(filename)) {
  saveRDS(myeloid.list, paste0("~/Jacinta/Nasal/Final_RDS_Files/",filename[[i]],".rds"))
}

