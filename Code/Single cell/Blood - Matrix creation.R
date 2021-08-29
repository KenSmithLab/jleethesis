##Extracting matrix from Seurat

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

setwd("~/Jacinta/Blood/Pre_processing/")

load("myeloid_final - Blood.RData")

#Create Seurat with all myeloid cells as "Myeloid"
new.cluster.ids <- c("Myeloid", "Myeloid", "Myeloid", "Myeloid")
names(new.cluster.ids) <- levels(myeloid.final)
myeloid.all <- RenameIdents(myeloid.final, new.cluster.ids)
table(myeloid.all@active.ident)


#Create Seurats with average gene expression for each sample
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
filename <- filename[grep("pbmc", filename)]
filename <- filename[-grep("cv15_pbmc_w12", filename)]
colnames(myeloid_all_matrix) <- filename

setwd("~/Jacinta/Blood/Module Preservation/")
save(myeloid_all_matrix, file = "myeloid_all_matrix.RData")


#Create Seurats with average gene expression for each sample and each myeloid cell type
myeloid_list <- SplitObject(myeloid.final, split.by = "sample")
myeloid_avg <- list()
for (i in 1:length(myeloid_list)) {
  myeloid_avg[[i]] <- AverageExpression(object = myeloid_list[[i]],
                    assays = "RNA", 
                    slot = "data", 
                    return.seurat = TRUE)
}



myeloid_avg[[1]]@assays$RNA@counts[1:30,]

myeloid_avg_merge <- merge(x=myeloid_avg[[1]],
                           y=c(myeloid_avg[[2]],myeloid_avg[[3]],myeloid_avg[[4]],myeloid_avg[[5]],myeloid_avg[[6]],
                               myeloid_avg[[7]],myeloid_avg[[8]],myeloid_avg[[9]],myeloid_avg[[10]],myeloid_avg[[11]],
                               myeloid_avg[[12]],myeloid_avg[[13]],myeloid_avg[[14]],myeloid_avg[[15]],myeloid_avg[[16]],
                               myeloid_avg[[17]],myeloid_avg[[18]],myeloid_avg[[19]],myeloid_avg[[20]],myeloid_avg[[21]],
                               myeloid_avg[[22]],myeloid_avg[[23]],myeloid_avg[[24]],myeloid_avg[[25]],myeloid_avg[[26]]))

myeloid_matrix <- GetAssayData(myeloid_avg_merge)

classical <- myeloid_matrix[,grep("Classical", colnames(myeloid_matrix))]
classical <- classical[,-grep("Non", colnames(classical))]
colnames(classical) <- filename

intermediate <- myeloid_matrix[,grep("Intermediate", colnames(myeloid_matrix))]
colnames(intermediate) <- filename

nonclassical <- myeloid_matrix[,grep("Non-Classical", colnames(myeloid_matrix))]
colnames(nonclassical) <- filename

myeloid_dc <- myeloid_matrix[,grep("Myeloid", colnames(myeloid_matrix))]
colnames(myeloid_dc) <- filename

save(classical, intermediate, nonclassical, myeloid_dc, file = "Myeloid_subsets_matrix.RData")


#Create matrices for all cells of each Myeloid subset ready for module preservation
classical_cells <- subset(myeloid.final, idents = c("Classical"))
classical_matrix <- 
  t(as.matrix(GetAssayData(classical_cells)))

intermediate_cells <- subset(myeloid.final, idents = c("Intermediate"))
intermediate_matrix <- 
  t(as.matrix(GetAssayData(intermediate_cells)))

nonclassical_cells <- subset(myeloid.final, idents = c("Non-Classical"))
nonclassical_matrix <- 
  t(as.matrix(GetAssayData(nonclassical_cells)))

myeloid_dc_cells <- subset(myeloid.final, idents = c("Myeloid DC"))
myeloid_dc_matrix <- 
  t(as.matrix(GetAssayData(myeloid_dc_cells)))

save(classical_cells, intermediate_cells, nonclassical_cells, myeloid_dc_cells,
     file = "Myeloid_subsets_seurat.RData")

save(classical_matrix, intermediate_matrix, nonclassical_matrix, myeloid_dc_matrix,
     file = "Myeloid_subsets_cells_matrix.RData")


##Subsetting classical monocytes and creating a matrix for each classical monocyte subpopulation
load("~/Jacinta/Blood/Pre-processing/myeloid_semifinal - Blood.RData")
DefaultAssay(myeloid.semifinal) <- "RNA"

new.cluster.ids <- c("zero", "one", "two", "three", "four", 
                     "five", "six", "seven", "eight", "nine", 
                     "ten", "eleven", "twelve", "thirteen")

names(new.cluster.ids) <- levels(myeloid.semifinal)
myeloid.semifinal <- RenameIdents(myeloid.semifinal, new.cluster.ids)

classical <- subset(myeloid.semifinal, idents = c("zero", "one", "two", "three", "four", 
                                                  "seven", "eight", "ten", "eleven", "thirteen"))
classical_list <- SplitObject(classical, split.by = "sample")

myeloid_avg <- list()
for (i in 1:length(classical_list)) {
  myeloid_avg[[i]] <- AverageExpression(object = classical_list[[i]],
                                        assays = "RNA", 
                                        slot = "data", 
                                        return.seurat = TRUE)
}

classical_avg_merge <- merge(x=myeloid_avg[[1]],
                           y=c(myeloid_avg[[2]],myeloid_avg[[3]],myeloid_avg[[4]],myeloid_avg[[5]],myeloid_avg[[6]],
                               myeloid_avg[[7]],myeloid_avg[[8]],myeloid_avg[[9]],myeloid_avg[[10]],myeloid_avg[[11]],
                               myeloid_avg[[12]],myeloid_avg[[13]],myeloid_avg[[14]],myeloid_avg[[15]],myeloid_avg[[16]],
                               myeloid_avg[[17]],myeloid_avg[[18]],myeloid_avg[[19]],myeloid_avg[[20]],myeloid_avg[[21]],
                               myeloid_avg[[22]],myeloid_avg[[23]],myeloid_avg[[24]],myeloid_avg[[25]],myeloid_avg[[26]]))

classical_matrix <- GetAssayData(classical_avg_merge)

classical_0 <- classical_matrix[, grep("zero", colnames(classical_matrix))]
colnames(classical_0) <- filename

classical_1 <- classical_matrix[, grep("one", colnames(classical_matrix))]
colnames(classical_1) <- filename

classical_2 <- classical_matrix[, grep("two", colnames(classical_matrix))]
colnames(classical_2) <- filename

classical_3 <- classical_matrix[, grep("three", colnames(classical_matrix))]
colnames(classical_3) <- filename

classical_4 <- classical_matrix[, grep("four", colnames(classical_matrix))]
colnames(classical_4) <- filename

classical_7 <- classical_matrix[, grep("seven", colnames(classical_matrix))]
colnames(classical_7) <- filename[c(1, 2, 6, 11, 15, 16, 18, 19, 24, 26)]

classical_8 <- classical_matrix[, grep("eight", colnames(classical_matrix))]
colnames(classical_8) <- filename[-21]

classical_10 <- classical_matrix[, grep("ten", colnames(classical_matrix))]
colnames(classical_10) <- filename[-18]

classical_11 <- classical_matrix[, grep("eleven", colnames(classical_matrix))]
colnames(classical_11) <- filename[c(4, 5, 6, 7, 9, 10, 12, 13,
                                     14, 15, 16, 18, 21, 22, 24)]

classical_13 <- classical_matrix[, grep("thirteen", colnames(classical_matrix))]
colnames(classical_13) <- filename[c(5, 9, 12, 14, 15, 17, 20, 26)]

samples_clusters <- as.data.frame(table(myeloid.semifinal$sample, myeloid.semifinal@active.ident))
samples_clusters <- samples_clusters %>% 
  pivot_wider(names_from = Var2,
              values_from = Freq)
write.csv(samples_clusters, "samples_myeloid_clusters.csv")

save(classical_0, classical_1, classical_2, classical_3, classical_4, 
     classical_7, classical_8, classical_10, classical_11, classical_13, 
     file = "Classical_subsets_avg_matrix.RData")


#Create Seurats with average gene expression for each sample of each timepoint then extract matrix
control <- myeloid_all_matrix[,grep("con", colnames(myeloid_all_matrix))]

zero <- myeloid_all_matrix[,grep("d1", colnames(myeloid_all_matrix))]

three <- myeloid_all_matrix[,grep("w12", colnames(myeloid_all_matrix))]

save(control, zero, three, file = "Myeloid_timepoints_matrix.RData")
