## Blood - Cluster annotations

library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(tidyverse)

setwd("~/Jacinta/Blood/Pre_processing/")

#Identifying broad clusters
load("VascSC_int_clustered - Blood.RData")
DefaultAssay(VascSC.int.clustered) <- "RNA"

pdf("Blood - Violin_Marker_plot1.pdf", width = 15, height = 13)
VlnPlot(VascSC.int.clustered, features = c("IL7R", "CCR7", "S100A4", "CD14", "LYZ", 
                                           "MS4A1", "CD8A", "FCGR3A", "MS4A7"))
dev.off()

pdf("Blood - Violin_Marker_plot2.pdf", width = 15, height = 13)
VlnPlot(VascSC.int.clustered, features = c("GNLY", "NKG7", "FCER1A", "CST3", "PPBP", 
                                           "S100A12", "SELL", "CD3E", "PRF1"))
dev.off()

pdf("Blood - Violin_Marker_plot3.pdf", width = 15, height = 13)
VlnPlot(VascSC.int.clustered, features = c("GPR183", "GZMA", "FOXP3", "IL2RA", 
                                           "SLC4A10", "TRAV1-2", "TRGV9", "TRDV2", 
                                           "KLRF1"))
dev.off()

pdf("Blood - Canonical_markers.pdf", width = 15, height = 13)
FeaturePlot(VascSC.int.clustered, features = c("CD3D", "CD3E", "MS4A1", "GNLY", "NKG7", "KLRF1",
                                               "CD14", "LYZ", "CD1C", "FCER1A", "PPBP", "HBB"), 
            order = TRUE)
dev.off()

#Renaming clusters
new.cluster.ids <- c("T", "T", "Myeloid", "T", "NK", "Myeloid", "B", "T", "Platelet", "HBB",
                     "Doublet", "T", "Myeloid", "NK", "B", "Myeloid", "Myeloid", "Myeloid", "Myeloid")

names(new.cluster.ids) <- levels(VascSC.int.clustered)
VascSC.int.clustered <- RenameIdents(VascSC.int.clustered, new.cluster.ids)

Idents(VascSC.int.clustered) <- factor(Idents(VascSC.int.clustered), 
                                  levels = c("Doublet", "HBB", "Platelet", "Myeloid", "NK", "B", "T"))
markers.to.plot <- c("CD3D", "CD3E", "MS4A1", "GNLY", "NKG7", "KLRF1",
                     "CD14", "LYZ", "PPBP", "HBB")
pdf("Blood - Canonical markers dotplot.pdf", height = 4, width = 6)
DotPlot(VascSC.int.clustered, features = markers.to.plot, dot.scale = 8) + 
  RotatedAxis()
dev.off()

pdf("Blood - Immune subsets.pdf", width = 5, height = 4)
DimPlot(VascSC.int.clustered, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

#No. of cells per cell type
table(VascSC.int.clustered@active.ident, VascSC.int.clustered$sample)

int.cell.time <- as.data.frame(table(VascSC.int.clustered@active.ident, VascSC.int.clustered$timepoint))

pdf(file = "Blood - Immune cell percent over time.pdf", width = 5, height = 3)
ggplot(int.cell.time, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  scale_x_discrete(labels = c("HC", "Relapse", "Remisson")) +
  theme_classic() +
  labs(x = element_blank(), y = "Portion of cells", fill = "Cell Type")
dev.off()





#Extract Myeloid Cells
myeloid.subset <- subset(x = VascSC.int.clustered, idents = "Myeloid")

save(myeloid.subset, file = "myeloid_subset.RData")

DefaultAssay(myeloid.subset) <- "RNA"

#Renormalise samples and integrate
myeloidList <- SplitObject(myeloid.subset, split.by = "sample")

myeloidList <- lapply(X = myeloidList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = myeloidList)
myeloidList <- lapply(X = myeloidList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

myeloid.anchors <- FindIntegrationAnchors(object.list = myeloidList,
                                         anchor.features = features,
                                         reduction = "rpca")
save(myeloid.anchors, file = "Myeloid_anchors - Blood.RData")

load("Myeloid_anchors - Blood.RData")
myeloid.integrated <- IntegrateData(anchorset = myeloid.anchors)
save(myeloid.integrated, file = "Myeloid_integrated - Blood.RData")

DefaultAssay(myeloid.integrated) <- "integrated"

myeloid.integrated <- ScaleData(myeloid.integrated, verbose = FALSE)
myeloid.integrated <- RunPCA(myeloid.integrated, npcs = 30, verbose = FALSE)

#Determine how many PCs to use, select around inflection point
ElbowPlot(myeloid.integrated)

myeloid.integrated <- RunUMAP(myeloid.integrated, reduction = "pca", dims = 1:17)
myeloid.integrated <- FindNeighbors(myeloid.integrated, reduction = "pca", dims = 1:17)
myeloid.clustered <- FindClusters(myeloid.integrated, resolution = 0.3)

save(myeloid.clustered, file = "Myeloid_clustered - Blood.RData")

#Re-Annotate
DefaultAssay(myeloid.clustered) <- "RNA"

pdf("Blood - Myeloid Canonical_markers.pdf", width = 15, height = 13)
FeaturePlot(myeloid.clustered, features = c("CD8A", "CD3E", "MS4A1", "GNLY", "NKG7", "KLRF1",
                                            "CD14", "LYZ", "CD1C", "FCER1A", "PPBP", "HBB"), 
            order = TRUE)
dev.off()

markers.to.plot2 <- c("CD8A", "CD3E", "MS4A1", "GNLY", "NKG7", "KLRF1",
                     "CD14", "CD1C", "FCER1A","LYZ", "PPBP", "HBB")
pdf("~/Jacinta/Blood/Pre_processing/Blood - Canonical markers dotplot2.pdf", height = 4, width = 6)
DotPlot(myeloid.clustered, features = markers.to.plot2, dot.scale = 8) + 
  RotatedAxis()
dev.off()

##Plot UMAP
pdf(file = "Blood - Initial Myeloid UMAP.pdf", width = 5, height = 4)
DimPlot(myeloid.clustered, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

pdf(file = "Blood - Initial Myeloid UMAP - Timepoint.pdf", width = 15, height = 5)
DimPlot(myeloid.clustered, reduction = "umap", split.by = "timepoint", raster = FALSE)
dev.off()

pdf(file = "Blood - Initial Myeloid UMAP - Batch.pdf", width = 15, height = 5)
DimPlot(myeloid.clustered, reduction = "umap", group.by = "batch", raster = FALSE)
dev.off()


#Subsetting to clear doublets, then re-normalising and integration
myeloid.trim <- subset(x = myeloid.clustered, idents = c(0, 1, 2, 3, 5, 6, 9, 10))

DefaultAssay(myeloid.trim) <- "RNA"

#Renormalise samples and integrate
myeloidtrimList <- SplitObject(myeloid.trim, split.by = "sample")

myeloidtrimList <- lapply(X = myeloidtrimList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = myeloidtrimList)
myeloidtrimList <- lapply(X = myeloidtrimList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

myeloid.trim.anchors <- FindIntegrationAnchors(object.list = myeloidtrimList,
                                          anchor.features = features, reduction = "rpca")
save(myeloid.trim.anchors, file = "Myeloid_trim_anchors - Blood.RData")

load("Myeloid_trim_anchors - Blood.RData")
myeloid.trim.integrated <- IntegrateData(anchorset = myeloid.trim.anchors, k.weight = 50)
save(myeloid.trim.integrated, file = "Myeloid_trim_integrated - Blood.RData")

DefaultAssay(myeloid.trim.integrated) <- "integrated"

myeloid.trim.integrated <- ScaleData(myeloid.trim.integrated, verbose = FALSE)
myeloid.trim.integrated <- RunPCA(myeloid.trim.integrated, npcs = 30, verbose = FALSE)

#Determine how many PCs to use, select around inflection point
ElbowPlot(myeloid.trim.integrated)

myeloid.trim.integrated <- RunUMAP(myeloid.trim.integrated, reduction = "pca", dims = 1:13)
myeloid.trim.integrated <- FindNeighbors(myeloid.trim.integrated, reduction = "pca", dims = 1:13)
save(myeloid.trim.integrated, file = "Myeloid_trim_int_before_cluster - Blood.RData")
myeloid.trim.clustered <- FindClusters(myeloid.trim.integrated, resolution = 0.7)

myeloid.trim.clustered2 <- FindClusters(myeloid.trim.integrated, resolution = 0.5, algorithm = 2)
myeloid.trim.clustered3 <- FindClusters(myeloid.trim.integrated, resolution = 0.5, algorithm = 3)

save(myeloid.trim.clustered, file = "Myeloid_trim_clustered - Blood.RData")

#Re-Annotate
load("Myeloid_trim_clustered - Blood.RData")
DefaultAssay(myeloid.trim.clustered) <- "RNA"

pdf("Blood - Myeloid_trim_Canonical_markers.pdf", width = 15, height = 13)
FeaturePlot(myeloid.trim.clustered, features = c("CD8A", "CD3E", "MS4A1", "GNLY", "NKG7", "KLRF1",
                                            "CD14", "LYZ", "CD1C", "FCER1A", "PPBP", "HBB"), 
            order = TRUE)
dev.off()

#Find cluster markers
all_markers <-FindAllMarkers(myeloid.trim.clustered, 
                             min.pct =  0.25, 
                             min.diff.pct = 0.25)
write.csv(all_markers, "Blood - Allmyeloidmarkers.csv")

top10markers <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
write.csv(top10markers, "Blood - Top10myeloidmarkers.csv")

##Plot UMAP
pdf(file = "Blood - Trimmed Myeloid UMAP.pdf", width = 10, height = 8)
DimPlot(myeloid.trim.clustered, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

pdf(file = "Blood - Trimmed Myeloid UMAP - Timepoint.pdf", width = 15, height = 5)
DimPlot(myeloid.trim.clustered, reduction = "umap", group.by = "timepoint", raster = FALSE)
dev.off()

pdf(file = "Blood - Trimmed Myeloid UMAP - Batch.pdf", width = 15, height = 5)
DimPlot(myeloid.trim.clustered, reduction = "umap", group.by = "batch", raster = FALSE)
dev.off()

#SingleR to tease out myeloid population
library(celldex)
library(SingleR)

monaco <- celldex::MonacoImmuneData()

monacoResults <- SingleR(test = as.SingleCellExperiment(myeloid.trim.clustered), ref = monaco, labels = monaco$label.fine)
write.csv(monacoResults, "monacoResults.csv")

results.test <- read.csv('monacoResults.csv', header = T)
rownames(results.test) = results.test$X
myeloid.trim.clustered <- AddMetaData(myeloid.trim.clustered, results.test)

DimPlot(myeloid.trim.clustered, group.by = 'pruned.labels', label = T, label.size = 5, repel = TRUE, raster = FALSE)

table(myeloid.trim.clustered$pruned.labels)

myeloid.semifinal <- myeloid.trim.clustered[,grep("monocytes|Myeloid", myeloid.trim.clustered$pruned.labels)]


pdf("UMAP split by SingleR labels.pdf", width = 15, height = 5)
DimPlot(myeloid.semifinal, split.by = 'pruned.labels', label = T, label.size = 5, repel = TRUE, raster = FALSE)
dev.off()

setwd("~/Jacinta/gene_sets/")
genelist <- list.files(full.names = TRUE, recursive = T, patt ="*.txt")
genelistfinal <- lapply(genelist, read.delim)
setwd("~/Jacinta")

CD16neg_vs_pos_UP <- genelistfinal[[1]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(CD16neg_vs_pos_UP), name = "CD16neg_vs_pos_UP")
M1_vs_M2_UP <- genelistfinal[[2]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(M1_vs_M2_UP), name = "M1_vs_M2_UP")
Mono_vs_M1_UP <- genelistfinal[[3]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(Mono_vs_M1_UP), name = "Mono_vs_M1_UP")
Mono_vs_M2_UP <- genelistfinal[[4]]
myeloid.semifinal <- AddModuleScore(myeloid.semifinal, list(Mono_vs_M2_UP), name = "Mono_vs_M2_UP")


p1 <- FeaturePlot(myeloid.semifinal, features = "CD16neg_vs_pos_UP1", order = TRUE)
p2 <- FeaturePlot(myeloid.semifinal, features = "M1_vs_M2_UP1", order = TRUE)
p3 <- FeaturePlot(myeloid.semifinal, features = "Mono_vs_M1_UP1", order = TRUE)
p4 <- FeaturePlot(myeloid.semifinal, features = "Mono_vs_M2_UP1", order = TRUE)
p5 <- p1|p2|p3|p4
p5

pdf("GeneSets in trimmed myeloid plot.pdf", width = 15, height = 5)
print(p5)
dev.off()

setwd("~/Jacinta/Blood/Pre_processing/")
pdf("Blood - Myeloid_markers.pdf", width = 15, height = 13)
FeaturePlot(myeloid.semifinal, features = c("CD14", "LYZ", "CD36",
                                                 "FCGR3A", "MS4A7","CCR5",
                                                 "CD86", "CD80", "FCGR1A", 
                                                 "SELPLG", "MRC1", "ARG1"), 
            order = TRUE)
dev.off()


save(myeloid.semifinal, file = "myeloid_semifinal - Blood.RData")

pdf(file = "Blood - Semifinal Myeloid UMAP.pdf", width = 10, height = 8)
DimPlot(myeloid.semifinal, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

new.cluster.ids <- c("Classical", "Classical", "Classical", "Classical", "Classical", 
                     "Non-Classical", "Intermediate", "Classical", "Classical", "Myeloid DC", 
                     "Classical", "Classical", "Non-Classical", "Classical")

names(new.cluster.ids) <- levels(myeloid.semifinal)
myeloid.final <- RenameIdents(myeloid.semifinal, new.cluster.ids)
table(myeloid.final$sample)

#Stacked bar plots for each timepoint
celltype.timepoint <- as.data.frame(table(myeloid.final@active.ident, myeloid.final$timepoint))

pdf(file = "Blood - Myeloid percent over time.pdf", width = 5, height = 3)
ggplot(celltype.timepoint, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  scale_x_discrete(labels = c("HC", "Relapse", "Remission")) +
  theme_classic() +
  labs(x = element_blank(), y = "Portion of cells", fill = "Cell Type")
dev.off()
  

all_markers <-FindAllMarkers(myeloid.semifinal, 
                             min.pct =  0.25, 
                             min.diff.pct = 0.25)
write.csv(all_markers, "All_markers - Blood.csv")

top10markers <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
write.csv(top10markers, "Top10markers - Blood.csv")

save(myeloid.final, file = "myeloid_final - Blood.RData")
saveRDS(myeloid.final, file = "myeloid_final.rds")

#Save individual RDS
demo <- read.csv("~/Jacinta/single_cell_phenodata.csv")
filename <- demo$filename
filename <- filename[grep("pbmc", filename)]
filename <- filename[-grep("cv15_pbmc_w12", filename)]

myeloid.list <- SplitObject(myeloid.final, split.by = "sample")
for (i in 1:length(filename)) {
  saveRDS(myeloid.list, paste0("~/Jacinta/Blood/Final_RDS_Files/",filename[[i]],".rds"))
}


#Plot final data
pdf(file = "Blood - Final Myeloid UMAP.pdf", width = 6, height = 4)
DimPlot(myeloid.final, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

pdf(file = "Blood - Final Myeloid UMAP - Timepoint.pdf", width = 15, height = 5)
DimPlot(myeloid.final, reduction = "umap", split.by = "timepoint", raster = FALSE)
dev.off()

pdf(file = "Blood - Final Myeloid UMAP - Batch.pdf", width = 15, height = 5)
DimPlot(myeloid.final, reduction = "umap", group.by = "batch", raster = FALSE)
dev.off()


#Plot data before renaming, splitting by timepoint
load("myeloid_semifinal - Blood.RData")
pdf(file = "Blood - Semifinal Myeloid UMAP - Timepoint.pdf", width = 15, height = 5)
DimPlot(myeloid.semifinal, reduction = "umap", split.by = "timepoint", raster = FALSE)
dev.off()

table(myeloid.semifinal@active.ident, myeloid.semifinal@meta.data$timepoint)

#Select classical subset and UMAP
classical_subsets <- subset(myeloid.semifinal, idents = c("0", "1", "2", "3", "4", "7", "8", "10", "11", "13"))
save(classical_subsets, file = "~/Jacinta/Blood/Classical Subset/classical_subsets.RData")

pdf(file = "Blood - Classical subsets UMAP.pdf", width = 5, height = 4)
DimPlot(classical_subsets, reduction = "umap", raster = FALSE, label = TRUE)
dev.off()

pdf(file = "Blood - Classical subsets UMAP - timepoint.pdf", width = 15, height = 5)
DimPlot(classical_subsets, reduction = "umap", split.by = "timepoint", 
        raster = FALSE, label = TRUE)
dev.off()

subsets_samples <- table(classical_subsets@active.ident, classical_subsets$sample)
table(classical_subsets$sample)

