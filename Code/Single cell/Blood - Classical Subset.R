##Classical subset analysis

library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(tidyverse)
library(fgsea)
library(data.table)

setwd("~/Jacinta/Blood/Classical Subset/")

load("classical_subsets.RData")

all_markers <- FindAllMarkers(classical_subsets)
write.csv(all_markers, "all_markers.csv")

#GSEA on each cluster based on findallmarkers
pathways.hallmark <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

all_markers <- read.csv("all_markers.csv")

#Cluster 0
cl0 <- all_markers[grep("0", all_markers$cluster),]
cl0 <- cl0[-grep("10", cl0$cluster),]
cl0 <- select(cl0, gene, avg_log2FC)
cl0 <- deframe(cl0)
cl0_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl0, nperm = 1000)

topPathwaysUp_custom <- cl0_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl0_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 0 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl0, cl0_fgsea, 
              gseaParam = 0.5)
dev.off()

cl0_fgsea_table <- cl0_fgsea[order(pval), ]
cl0_fgsea_table <- as.data.frame(cl0_fgsea_table)
cl0_fgsea_table <- apply(cl0_fgsea_table,2,as.character)
write.csv(cl0_fgsea_table, file = "cl0_fgsea_table.csv")

plotEnrichment(pathways.hallmark[["HALLMARK_COMPLEMENT"]],
               cl0) + 
    labs(title = "HALLMARK_COMPLEMENT")


#Cluster 1
cl1 <- all_markers[grep("1", all_markers$cluster),]
cl1 <- cl1[-grep(c("10|11|13"), cl1$cluster),]
cl1 <- select(cl1, gene, avg_log2FC)
cl1 <- deframe(cl1)
cl1_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl1, nperm = 1000)

topPathwaysUp_custom <- cl1_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl1_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 1 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl1, cl1_fgsea, 
              gseaParam = 0.5)
dev.off()

cl1_fgsea_table <- cl1_fgsea[order(pval), ]
cl1_fgsea_table <- as.data.frame(cl1_fgsea_table)
cl1_fgsea_table <- apply(cl1_fgsea_table,2,as.character)
write.csv(cl1_fgsea_table, file = "cl1_fgsea_table.csv")

#Cluster 2
cl2 <- all_markers[grep("2", all_markers$cluster),]
cl2 <- select(cl2, gene, avg_log2FC)
cl2 <- deframe(cl2)
cl2_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl2, nperm = 1000)

topPathwaysUp_custom <- cl2_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl2_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 2 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl2, cl2_fgsea, 
              gseaParam = 0.5)
dev.off()

cl2_fgsea_table <- cl2_fgsea[order(pval), ]
cl2_fgsea_table <- as.data.frame(cl2_fgsea_table)
cl2_fgsea_table <- apply(cl2_fgsea_table,2,as.character)
write.csv(cl2_fgsea_table, file = "cl2_fgsea_table.csv")

#Cluster 3
cl3 <- all_markers[grep("3", all_markers$cluster),]
cl3 <- cl3[-grep("13", cl3$cluster),]
cl3 <- select(cl3, gene, avg_log2FC)
cl3 <- deframe(cl3)
cl3_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl3, nperm = 1000)

topPathwaysUp_custom <- cl3_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl3_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 3 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl3, cl3_fgsea, 
              gseaParam = 0.5)
dev.off()

cl3_fgsea_table <- cl3_fgsea[order(pval), ]
cl3_fgsea_table <- as.data.frame(cl3_fgsea_table)
cl3_fgsea_table <- apply(cl3_fgsea_table,2,as.character)
write.csv(cl3_fgsea_table, file = "cl3_fgsea_table.csv")

#Cluster 4
cl4 <- all_markers[grep("4", all_markers$cluster),]
cl4 <- select(cl4, gene, avg_log2FC)
cl4 <- deframe(cl4)
cl4_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl4, nperm = 1000)

topPathwaysUp_custom <- cl4_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl4_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 4 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl4, cl4_fgsea, 
              gseaParam = 0.5)
dev.off()

cl4_fgsea_table <- cl4_fgsea[order(pval), ]
cl4_fgsea_table <- as.data.frame(cl4_fgsea_table)
cl4_fgsea_table <- apply(cl4_fgsea_table,2,as.character)
write.csv(cl4_fgsea_table, file = "cl4_fgsea_table.csv")

#Cluster 7
cl7 <- all_markers[grep("7", all_markers$cluster),]
cl7 <- select(cl7, gene, avg_log2FC)
cl7 <- deframe(cl7)
cl7_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl7, nperm = 1000)

topPathwaysUp_custom <- cl7_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl7_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 7 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl7, cl7_fgsea, 
              gseaParam = 0.5)
dev.off()

cl7_fgsea_table <- cl7_fgsea[order(pval), ]
cl7_fgsea_table <- as.data.frame(cl7_fgsea_table)
cl7_fgsea_table <- apply(cl7_fgsea_table,2,as.character)
write.csv(cl7_fgsea_table, file = "cl7_fgsea_table.csv")

#Cluster 8
cl8 <- all_markers[grep("8", all_markers$cluster),]
cl8 <- select(cl8, gene, avg_log2FC)
cl8 <- deframe(cl8)
cl8_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl8, nperm = 1000)

topPathwaysUp_custom <- cl8_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl8_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 8 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl8, cl8_fgsea, 
              gseaParam = 0.5)
dev.off()

cl8_fgsea_table <- cl8_fgsea[order(pval), ]
cl8_fgsea_table <- as.data.frame(cl8_fgsea_table)
cl8_fgsea_table <- apply(cl8_fgsea_table,2,as.character)
write.csv(cl8_fgsea_table, file = "cl8_fgsea_table.csv")

#Cluster 10
cl10 <- all_markers[grep("10", all_markers$cluster),]
cl10 <- select(cl10, gene, avg_log2FC)
cl10 <- deframe(cl10)
cl10_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl10, nperm = 1000)

topPathwaysUp_custom <- cl10_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl10_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 10 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl10, cl10_fgsea, 
              gseaParam = 0.5)
dev.off()

cl10_fgsea_table <- cl10_fgsea[order(pval), ]
cl10_fgsea_table <- as.data.frame(cl10_fgsea_table)
cl10_fgsea_table <- apply(cl10_fgsea_table,2,as.character)
write.csv(cl10_fgsea_table, file = "cl10_fgsea_table.csv")

#Cluster 11
cl11 <- all_markers[grep("11", all_markers$cluster),]
cl11 <- select(cl11, gene, avg_log2FC)
cl11 <- deframe(cl11)
cl11_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl11, nperm = 1000)

topPathwaysUp_custom <- cl11_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl11_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 11 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl11, cl11_fgsea, 
              gseaParam = 0.5)
dev.off()

cl11_fgsea_table <- cl11_fgsea[order(pval), ]
cl11_fgsea_table <- as.data.frame(cl11_fgsea_table)
cl11_fgsea_table <- apply(cl11_fgsea_table,2,as.character)
write.csv(cl11_fgsea_table, file = "cl11_fgsea_table.csv")

#Cluster 13
cl13 <- all_markers[grep("13", all_markers$cluster),]
cl13 <- select(cl13, gene, avg_log2FC)
cl13 <- deframe(cl13)
cl13_fgsea <- fgsea(pathways = pathways.hallmark, stats = cl13, nperm = 1000)

topPathwaysUp_custom <- cl13_fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- cl13_fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Cluster 13 - GSEA table.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_custom], cl13, cl13_fgsea, 
              gseaParam = 0.5)
dev.off()

cl13_fgsea_table <- cl13_fgsea[order(pval), ]
cl13_fgsea_table <- as.data.frame(cl13_fgsea_table)
cl13_fgsea_table <- apply(cl13_fgsea_table,2,as.character)
write.csv(cl13_fgsea_table, file = "cl13_fgsea_table.csv")


##EnrichR for each cluster
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

#Positive genes
#Create list of positive genes for each cluster
pos_markers <- all_markers %>% 
    subset(avg_log2FC > 0)
pos_markers <- split(pos_markers, pos_markers$cluster)

clusters <- c("0", "1", "2", "3", "4", "7", "8", "10", "11", "13")
gene_vector <- list()
enriched_pos <- list()
hallmark_pos <- list()
for (i in clusters) {
    gene_vector[[i]] <- as.vector(pos_markers[[i]]$gene)
    enriched_pos[[i]] <- enrichr(gene_vector[[i]], "MSigDB_Hallmark_2020")
    hallmark_pos[[i]] <- enriched_pos[[i]][["MSigDB_Hallmark_2020"]]
    write.csv(hallmark_pos[[i]], paste0("Enrichr - Positive Hallmark - Cluster", i, ".csv"))
}

#Negative genes
#Create list of positive genes for each cluster
neg_markers <- all_markers %>% 
    subset(avg_log2FC < 0)
neg_markers <- split(neg_markers, neg_markers$cluster)

gene_vector <- list()
enriched_neg <- list()
hallmark_neg <- list()
for (i in clusters) {
    gene_vector[[i]] <- as.vector(neg_markers[[i]]$gene)
    enriched_neg[[i]] <- enrichr(gene_vector[[i]], "MSigDB_Hallmark_2020")
    hallmark_neg[[i]] <- enriched_neg[[i]][["MSigDB_Hallmark_2020"]]
    write.csv(hallmark_neg[[i]], paste0("Enrichr - Negative Hallmark - Cluster", i, ".csv"))
}



##Magenta Eigengene for each cluster
magentaGenes <- read.csv("~/Jacinta/magentaGenes.csv")
magentaGenes <- magentaGenes$x
pinkGenes <- read.csv("~/Jacinta/pinkGenes.csv")
pinkGenes <- pinkGenes$x

load("~/Jacinta/Blood/Pre_processing/Classical_subsets_avg_matrix.RData")
magenta <- rep("magenta", times = length(magentaGenes))
pink <- rep("pink", times = length(pinkGenes))

classical_list <- list(classical_0, classical_1, classical_2, classical_3,
                       classical_4, classical_7, classical_8, classical_10,
                       classical_11, classical_13)
clusters <- c("cl0", "cl1", "cl2", "cl3", "cl4", "cl7", "cl8", "cl10", "cl11", "cl13")
magenta_list <- list()
pink_list <- list()
magentaME_list <- list()
MEmagenta_list <- list()
pinkME_list <- list()
MEpink_list <- list()
magenta_name <- list()
pink_name <- list()
for (i in 1:length(classical_list)) {
    classical_list[[i]] <- as.data.frame(classical_list[[i]])
    magenta_list[[i]] <- classical_list[[i]][magentaGenes,]
    magenta_list[[i]] <- t(magenta_list[[i]])
    magentaME_list[[i]] <- moduleEigengenes(magenta_list[[i]], magenta)
    MEmagenta_list[[i]] <- magentaME_list[[i]]$eigengenes$MEmagenta
    magenta_name[[i]] <- rep(paste0(clusters[[i]]), times = length(MEmagenta_list[[i]]))
    MEmagenta_list[[i]] <- cbind(MEmagenta_list[[i]],colnames(classical_list[[i]]), magenta_name[[i]])
    
    pink_list[[i]] <- classical_list[[i]][pinkGenes,]
    pink_list[[i]] <- t(pink_list[[i]])
    pinkME_list[[i]] <- moduleEigengenes(pink_list[[i]], pink)
    MEpink_list[[i]] <- pinkME_list[[i]]$eigengenes$MEpink
    pink_name[[i]] <- rep(paste0(clusters[[i]]), times = length(MEpink_list[[i]]))
    MEpink_list[[i]] <- cbind(MEpink_list[[i]], colnames(classical_list[[i]]), pink_name[[i]])
}


##MAGENTA
mag_list <- map(MEmagenta_list, as.data.table)
mag <- rbindlist(mag_list, fill = TRUE)
colnames(mag) <- c("magentaME", "filename", "cluster")
mag$filename <- str_remove(mag$filename, "_cite")
mag_final <- mag %>% 
    separate(
        col=filename, into=c("id","tissue", "timepoint"), convert = TRUE, sep="_")
mag_final <- transform(mag_final, magentaME = as.numeric(magentaME))
mag_final$timepoint <- as.factor(mag_final$timepoint)
mag_final$cluster <- factor(mag_final$cluster, 
                             levels = c("cl0", "cl1", "cl2", "cl3", "cl4",
                                        "cl7", "cl8", "cl10", "cl11", "cl13"))

#Stats between pinkME and Time per cluster
mag_aov <- aov(magentaME ~ timepoint + cluster, data = mag_final)
# Summary of the analysis
summary(mag_aov)
#Determine which groups are significant
mag_bonf <- mag_final %>% 
    group_by(cluster) %>% 
    pairwise_t_test(
        magentaME ~ timepoint, pool.sd = FALSE
    )
mag_bonf <-mag_bonf %>% 
    add_xy_position(x = "timepoint")
mag_bonf

pdf(file = "Blood - Magenta module over time per classical cluster.pdf", height = 3, width = 9)
ggplot(mag_final, aes(x = cluster, y = magentaME, fill = timepoint)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("skyblue1", "brown1", "darkorange1"), name = "Timepoint",
                      labels = c("HC", "AAV Zero", "AAV Three")) +
    scale_x_discrete(labels = c("0", "1", "2", "3", "4", "7", "8", "10", "11", "13")) +
    theme_classic() +
    labs(x = "Cluster", y = "Magenta Eigengene") +
    ylim(-0.005, 0.008)
#stat_pvalue_manual(pink_bonf, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)
dev.off()


##PINK
pin_list <- map(MEpink_list, as.data.table)
pink <- rbindlist(pin_list, fill = TRUE)
colnames(pink) <- c("pinkME", "filename", "cluster")
pink$filename <- str_remove(pink$filename, "_cite")
pink_final <- pink %>% 
    separate(
        col=filename, into=c("id","tissue", "timepoint"), convert = TRUE, sep="_")
pink_final <- transform(pink_final, pinkME = as.numeric(pinkME))
pink_final$timepoint <- as.factor(pink_final$timepoint)
pink_final$cluster <- factor(pink_final$cluster, 
                             levels = c("cl0", "cl1", "cl2", "cl3", "cl4",
                                        "cl7", "cl8", "cl10", "cl11", "cl13"))

#Stats between pinkME and Time per cluster
pink_aov <- aov(pinkME ~ timepoint + cluster, data = pink_final)
# Summary of the analysis
summary(pink_aov)
#Determine which groups are significant
pink_bonf <- pink_final %>% 
    group_by(cluster) %>% 
    pairwise_t_test(
        pinkME ~ timepoint, pool.sd = FALSE
    )
pink_bonf <-pink_bonf %>% 
    add_xy_position(x = "timepoint")
pink_bonf

pdf(file = "Blood - Pink module over time per classical cluster.pdf", height = 3, width = 9)
ggplot(pink_final, aes(x = cluster, y = pinkME, fill = timepoint)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("skyblue1", "brown1", "darkorange1"), name = "Timepoint",
                      labels = c("HC", "AAV Zero", "AAV Three")) +
    scale_x_discrete(labels = c("0", "1", "2", "3", "4", "7", "8", "10", "11", "13")) +
    theme_classic() +
    labs(x = "Cluster", y = "Pink Eigengene") +
    ylim(-0.005,0.01)
    #stat_pvalue_manual(pink_bonf, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)
dev.off()

