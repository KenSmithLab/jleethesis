##Nasal - Extracting matrix from Seurat

library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(DESeq2)

setwd("~/Jacinta/Nasal/Pre_processing/")

load("myeloid_semifinal - Nasal.RData")

#Create Seurat with all myeloid cells as "Myeloid"
new.cluster.ids <- c("Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid")
names(new.cluster.ids) <- levels(myeloid.semifinal)
myeloid.all <- RenameIdents(myeloid.semifinal, new.cluster.ids)
table(myeloid.all@active.ident)


#Create Seurats with average gene expression for each sample then extract matrix
myeloid_all_list <- SplitObject(myeloid.all, split.by = "sample")
myeloid_avg <- list()
matrix <- list()
for (i in 1:length(myeloid_all_list)) {
  myeloid_avg[[i]] <- AverageExpression(object = myeloid_all_list[[i]],
                                        assays = "RNA", 
                                        slot = "data", 
                                        return.seurat = TRUE)
  matrix[[i]] <- as.matrix(GetAssayData(myeloid_avg[[i]]))
}
myeloid_all_matrix <- as.data.frame(matrix)

demo <- read.csv("~/Jacinta/single_cell_phenodata.csv")
filename <- demo$filename
filename <- filename[grep("nasal", filename)]
filename <- filename[-c(3, 4, 6, 8, 9, 12, 19, 21, 23)]
colnames(myeloid_all_matrix) <- filename

setwd("~/Jacinta/Nasal/Module Preservation/")
save(myeloid_all_matrix, file = "myeloid_all_matrix.RData")

#Create Seurats with average gene expression for each sample of each timepoint then extract matrix
control <- myeloid_all_matrix[,grep("con", colnames(myeloid_all_matrix))]

zero <- myeloid_all_matrix[,grep("d1", colnames(myeloid_all_matrix))]

three <- myeloid_all_matrix[,grep("w12", colnames(myeloid_all_matrix))]

save(control, zero, three, file = "Myeloid_timepoints_matrix.RData")

