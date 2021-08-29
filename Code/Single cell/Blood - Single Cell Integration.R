##Single cell - Blood Integration

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(BiocManager)
library(metap)
library(multtest)
library(reticulate)

setwd("~/Jacinta/Blood/Initial_RDS_Files/")

files <- list.files()
files <- files[grep('pbmc', files)]
#Remove rds to get final sample name
sample.name <- substr(files, start=1,stop=nchar(files)-4)

#read in demographics

setwd("~/Jacinta/")
demo <- read.csv('~/Jacinta/single_cell_phenodata.csv', header = T)
demo <- demo[grep('pbmc', demo$file, invert = F),]

#File names matching with demo
file.names <- character()
for(i in 1:length(files)){
  file.names[i] <- files[grep(demo$file[i],files)]
}

#Read in RDS and add meta data
setwd("~/Jacinta/Blood/Initial_RDS_Files/")
counts <- list()
data <- list()
for (i in 1:length(file.names)) {
  data[[i]] <- readRDS(files[[i]])
  data[[i]]@meta.data$sample <- demo$filename[i]
  data[[i]]@meta.data$id <- demo$Trial_ID[i]
  data[[i]]@meta.data$timepoint <- demo$timepoint[i] 
  data[[i]]@meta.data$batch <- demo$batch[i]
  data[[i]]@meta.data$age <- demo$age[i]
  data[[i]]@meta.data$sex <- demo$sex[i]
  data[[i]]@meta.data$bvas <- demo$bvas_score[i]
  data[[i]]@meta.data$pr3 <- demo$pr3_titre[i]
  data[[i]]@meta.data$crp <- demo$crp[i]
  data[[i]]@meta.data$gfr <- demo$gfr[i]
  data[[i]]@meta.data$s_aureus <- demo$S_aureus[i]
  data[[i]]@meta.data$past_aza <- demo$past_aza[i]
  data[[i]]@meta.data$past_rtx <- demo$past_rtx[i]
  data[[i]]@meta.data$past_mtx <- demo$past_mtx[i]
  data[[i]]@meta.data$past_cyc <- demo$past_cyc[i]
  data[[i]]@meta.data$ent <- demo$ent[i]
  data[[i]]@meta.data$lung <- demo$lung[i]
  data[[i]]@meta.data$joints <- demo$joints[i]
  data[[i]]@meta.data$skin <- demo$skin[i]
  data[[i]]@meta.data$renal<- demo$renal[i]
  data[[i]]@meta.data$muscle <- demo$muscle[i]
  counts[[i]] <- length(colnames(data[[i]]))
  names(counts)[[i]] <- file.names[i]
}

setwd("~/Jacinta/Blood/Pre_processing/")
cellcounts <- do.call(rbind, counts)
x <- as.data.frame(cellcounts)
sum(x$V1)
write.csv(cellcounts,'Blood_counts.csv')

#Create Seurat objects for each file
for (i in 1:length(data)) {
  
  normcounts <- GetAssayData(data[[i]], slot = 'data')
  counts <- GetAssayData(data[[i]], slot = 'counts')
  metadata <-  data[[i]]@meta.data
  
  colnames(normcounts) <- paste0(data[[i]]$sample,"_", gsub('-.*','',colnames(normcounts)))
  colnames(counts) <- paste0(data[[i]]$sample,"_", gsub('-.*','',colnames(counts)))
  rownames(metadata) <- paste0(data[[i]]$sample,"_", gsub('-.*','',rownames(metadata)))
  
  
  data[[i]] <- CreateSeuratObject(
    counts,
    project = "AAV",
    assay = "RNA",
    min.cells = ,
    min.features = 0,
    names.field = 1,
    names.delim = "_",
    meta.data = metadata
  )
  
  data[[i]][["percent.mt"]] <- PercentageFeatureSet(data[[i]], pattern = "^MT-")
  data[[i]] <- subset(data[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
}

#Merge data
mergedData <- merge(x=data[[1]],
                    y=c(data[[2]], data[[3]], data[[4]],data[[5]], 
                        data[[6]], data[[7]], data[[8]], data[[9]],
                        data[[10]], data[[11]], data[[12]], data[[13]],
                        data[[14]], data[[15]], data[[16]], data[[17]],
                        data[[19]], data[[20]], data[[21]],
                        data[[22]], data[[23]], data[[24]],
                        data[[25]], data[[26]], data[[27]]))

save(mergedData, file = "Blood - MergedVascSC.RData")

load("~/Jacinta/Blood/Pre_processing/Blood - MergedVascSC.RData")
length(Cells(mergedData))

pdf("~/Jacinta/Blood/Pre_processing/Blood - QC violin plot.pdf", width = 15, height = 5)
VlnPlot(mergedData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

VascSClist <- SplitObject(mergedData, split.by = "sample")

#Use SCtransform to normalize data
#glmGamPoi method used to enable faster estimation of regression parameters in SCTransform
#VascSClist <- lapply(X = VascSClist, FUN = SCTransform, method = "glmGamPoi")

#Integrate Data using fast PCA method
#FPCA method is better with data for large number datasets
#features <- SelectIntegrationFeatures(object.list = VascSClist, nfeatures = 3000)
#VascSClist <- PrepSCTIntegration(object.list = VascSClist, anchor.features = features)
#VascSClist <- lapply(X = VascSClist, FUN = RunPCA, features = features)

#VascSC.anchors <- FindIntegrationAnchors(object.list = VascSClist, 
#                                         normalization.method = "SCT",
#                                         anchor.features = features,
#                                         dims = 1:30,
#                                         reduction = "rpca",
#                                         k.anchor = 20)
#save(VascSC.anchors, file = "VascSC_anchors - Blood.RData")

#VascSC.integrated <- IntegrateData(anchorset = VascSC.anchors, 
#                                   normalization.method = "SCT",
#                                   dims = 1:30)

#DefaultAssay(VascSC.integrated) <- "integrated"
#VascSC.integrated <- RunPCA(VascSC.integrated, verbose = FALSE)
#VascSC.integrated <- RunUMAP(VascSC.integrated, reduction = "pca", dims = 1:30)

#VascSC.integrated <- FindNeighbors(VascSC.integrated, dims = 1:30, verbose = FALSE)
#VascSC.integrated <- FindClusters(VascSC.integrated, verbose = FALSE)

#save(VascSC.integrated, file = "VascSCintegrated data - Blood.RData")


## Using normal way to normalize data
# normalize and identify variable features for each dataset independently
VascSClist <- lapply(X = VascSClist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = VascSClist)
VascSClist <- lapply(X = VascSClist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Integration
VascSC.anchors <- FindIntegrationAnchors(object.list = VascSClist,
                                         anchor.features = features,
                                         reduction = "rpca")
save(VascSC.anchors, file = "VascSC_anchors - Blood.RData")

load("VascSC_anchors - Blood.RData")
VascSC.integrated <- IntegrateData(anchorset = VascSC.anchors)
save(VascSC.integrated, file = "VascSC_integrated - Blood.RData")

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(VascSC.integrated) <- "integrated"

#Determine dimensionality of data
VascSC.integrated <- ScaleData(VascSC.integrated, verbose = FALSE)

VascSC.integrated <- RunPCA(VascSC.integrated, npcs = 30, verbose = FALSE)

pdf("~/Jacinta/Blood/Pre_processing/Blood - Elbowplot1.pdf", width = 5, height = 5)
ElbowPlot(VascSC.integrated)
dev.off()


# Run the standard workflow for visualization and clustering
VascSC.integrated <- RunUMAP(VascSC.integrated, reduction = "pca", dims = 1:19)
VascSC.integrated <- FindNeighbors(VascSC.integrated, reduction = "pca", dims = 1:19)
VascSC.int.clustered <- FindClusters(VascSC.integrated, resolution = 0.3)

save(VascSC.int.clustered, file = "VascSC_int_clustered - Blood.RData")

##Plot UMAP
pdf(file = "Blood - Initial UMAP plot after integration.pdf", width = 5, height = 4)
DimPlot(VascSC.int.clustered, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

pdf(file = "Blood - Initial UMAP plot after integration - Timepoint.pdf", width = 10, height = 8)
DimPlot(VascSC.int.clustered, reduction = "umap", group.by = "timepoint", raster = FALSE)
dev.off()

pdf(file = "Blood - Initial UMAP plot after integration - Batch.pdf", width = 5, height = 4)
DimPlot(VascSC.int.clustered, reduction = "umap", group.by = "batch", raster = FALSE) +
  labs(title = element_blank())
dev.off()

