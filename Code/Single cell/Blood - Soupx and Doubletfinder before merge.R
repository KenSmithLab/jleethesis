##Single cell data input using Soupx And doublet finder before merge
#Blood

library(Seurat)
library(patchwork)
library(DoubletFinder)
library(SoupX)
library(Matrix)
library(DropletUtils)
library(methods)

setwd("~/Jacinta/Blood/Blood_Raw_Data/")

samples <- list.files()
length(samples)

for(i in 1:length(samples)){
  
  toc = Read10X(paste0(samples[[i]],'/filtered_feature_bc_matrix/'))
  tod = Read10X(paste0(samples[[i]],'/raw_feature_bc_matrix/'))
  #sc = SoupChannel(tod$`Gene Expression`, toc$`Gene Expression`) # if CITE-seq
  sc = SoupChannel(tod, toc) # if no CITE-seq
  
  #seu <- CreateSeuratObject(counts = toc$`Gene Expression`) # if CITE-seq
  seu <- CreateSeuratObject(counts = toc) # if no CITE-seq
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 500)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- FindNeighbors(seu, k.param = 10, dims = 1:20)
  seu <- FindClusters(seu, algorithm = 1)
  seu <- RunUMAP(seu, dims = 1:20, n.neighbors = 10, min.dist = 0.3)
  clusters <- Idents(seu)
  
  seu_metaData = cbind(Embeddings(seu, reduction = "umap"), seu@meta.data)
  sc = setClusters(sc, setNames(seu_metaData$seurat_clusters, rownames(seu_metaData)))
  sc = setDR(sc, seu_metaData[colnames(seu), c("UMAP_1", "UMAP_2")])
  sc = autoEstCont(sc)
  out = adjustCounts(sc)
  DropletUtils:::write10xCounts(paste0(samples[[i]],'/strainedCounts'), out)
  
  
  #DoubletFinder
  data <- Read10X(data.dir = paste0(samples[[i]],'/strainedCounts'))
  data <- CreateSeuratObject(data, min.cells = 0, min.features = 200)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  #data <- subset(data, subset = percent.mt < 20 & nFeature_RNA < 2500)
  data <- NormalizeData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 500, verbose = FALSE)
  data <- ScaleData(data, vars.to.regress = c('percent.mt', 'nCount_RNA'), verbose = FALSE)
  data <- RunPCA(data, features = VariableFeatures(object = data))
  data <- RunUMAP(data, dims = 1:20, n.neighbors = 10, min.dist = 0.3)
  data <- FindNeighbors(data, k.param = 10, dims = 1:20)
  data <- FindClusters(data, algorithm = 1)
  
  #remove doublets
  
  sweep.res.list <- paramSweep_v3(data, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           
  nExp_poi <- round(0.075*length(names(data$orig.ident))) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  data <- doubletFinder_v3(data, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  pann <- colnames(data@meta.data)[grep("pANN",colnames(data@meta.data))]
  data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pann, sct = FALSE)
  df <- colnames(data@meta.data)[grep("DF",colnames(data@meta.data))]
  df <- df[grep(substr(pann,start=nchar(pann)-3,stop=nchar(pann)),df)]
  #DimPlot(data, reduction = "umap", group.by = df)
  col <- which(colnames(data@meta.data)==df)
  row <- grep("Singlet",data@meta.data[,col])
  data <- data[,row]
  
  saveRDS(data, paste0('~/Jacinta/Blood/Initial_RDS_Files/',samples[[i]],'.rds'))
}
