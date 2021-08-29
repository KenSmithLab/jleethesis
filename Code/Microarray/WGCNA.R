#WGCNA
#Weighted gene co-expression network analysis

library(oligo)
library(limma)
library(WGCNA)
library(affycoretools)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(xlsx)
library(RColorBrewer)

setwd("~/VascMicroarray/Usable Data")

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#Read in files that have been QCed for outliers
load(file = "Vascmicroarray-final.RData")

#Transpose
cleanMatrix <- t(vascFinalMatrix)

#Loading in numeric phenodata (traitdata)
traitdata <- read.csv("Phenodata_numeric.csv") 
dim(traitdata)
names(traitdata)
#traitdata <- traitdata[order(traitdata$filename), ]

cleantraits <- traitdata[, -1]

#Similarity matrix - Similarity between genes
#similarity = cor(cleanMatrix, method = "spearman")

#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(cleanMatrix, verbose = 3);
gsg$allOK
#Outcome is TRUE, so no need to remove outliers

# Cluster samples
wgcnaTree = hclust(dist(cleanMatrix), method = "average")

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(wgcnaTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


#Outliers can be removed in R, but removed in CEL file in original and preprocessing run again
# Plot a line to show the cut
#abline(h = 100, col = "red");
# Determine cluster under the line
#clust = cutreeStatic(wgcnaTree, cutHeight = 100, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#cleanMatrix = cleanMatrix[keepSamples, ]
nGenes = ncol(cleanMatrix)
nSamples = nrow(cleanMatrix)

#Remove outlier in traitdata, but removed in 
#traitdata <- traitdata[-160, ]
#cleantraits <- traitdata[, -1]

#Repeat clustering
wgcnaTree = hclust(dist(cleanMatrix), method = "average")
pdf(file = "Sample clustering.pdf", width = 16, height = 5)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(wgcnaTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# Convert traits to a color representation: white means low, red means high, 
#grey means missing entry
traitColours = numbers2colors(cleantraits, signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(wgcnaTree, traitColours,
                    groupLabels = names(cleantraits),
                    main = "Sample dendrogram and trait heatmap")

#save
save(cleanMatrix, cleantraits, finalAnnots, file = "VascWGCNA_DataInput.RData")


##Network Construction

# Define the adjacency matrix using soft thresholding with beta=6
#ADJ1=abs(cor(cleanMatrix,use="p"))^6
# When you have a lot of genes use the following code
k=softConnectivity(datE=cleanMatrix,power=6)
# Plot a histogram of k and a scale free topology plot
pdf(file = "Hist of k and Check scale free topology.pdf", width = 15, height = 8)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(cleanMatrix, powerVector = powers, verbose = 5)

sft$powerEstimate

pdf(file = "Scale independence - Mean connectivity.pdf", width = 16, height = 8)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#Power 9 used

#One-step network construction and module detection
net = blockwiseModules(cleanMatrix, power = 9,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Vascmodules",
                       verbose = 3)
table(net$colors)

# open a graphics window
pdf(file = "Cluster Dendrogram.pdf", width = 18, height = 10)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#Save
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, net,
     file = "VascNetworkConstruction.RData")


##Relating modules to clinical traits
load(file = "VascNetworkConstruction.RData")
# Define numbers of genes and samples
nGenes = ncol(cleanMatrix);
nSamples = nrow(cleanMatrix);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(cleanMatrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, cleantraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#Displaying in table
# Will display correlations and their p-values
pdf(file = "Module-trait relationships.pdf", width = 16, height = 15)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
VascHeatmap = labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(cleantraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#tAnnots = t(finalAnnots)

probes = colnames(cleanMatrix)
probes2annot = match(probes, finalAnnots$PROBEID)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

#Interfacing network analysis with other data
allEntrezID = finalAnnots$ENTREZID[probes2annot]

#Choose interesting modules in Vasculitis
#vascModules = c("midnightblue", "magenta", "salmon")
#for (module in vascModules)
#{
        # Select module probes
#        modGenes = (moduleColors==module)
        # Get their entrez ID codes
#        modEnIDs = allEntrezID[modGenes];
        # Write them into a file
#        fileName = paste("VascEntrezIDs-", module, ".txt", sep="");
#        write.table(as.data.frame(modEnIDs), file = fileName,
#                    row.names = FALSE, col.names = FALSE)
#}

# As background in the enrichment analysis, we will use all probes in the analysis.
#fileName = paste("EntrezID-all.txt", sep="");
#write.table(as.data.frame(allEntrezID), file = fileName,
#            row.names = FALSE, col.names = FALSE)

#salmonEntrez = read.table("VascEntrezIDs-salmon.txt")
#magentaEntrez = read.table("VascEntrezIDs-magenta.txt")
#midnightblueEntrez = scan("VascEntrezIDs-midnightblue.txt", 
#                          what = list(NULL, name = character()))


#Enrichment analysis
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)


GOenr = GOenrichmentAnalysis(moduleColors, allEntrezID, organism = "human", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

#Magenta, Pink for looking at between timepoints 
#Midnight blue for PR3 and MPO titre


#Network Visualization

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(cleanMatrix, power = 9);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
#diag(plotTOM) = NA;
# Call the plot function
#pdf(file = "Network heatmap plot - all genes.pdf", width = 15, height = 15)
#sizeGrWindow(9,9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#dev.off()

#Visualising network of eigengenes comparing active_vs_hc
active_vs_hc = as.data.frame(cleantraits$active_vs_hc)
names(active_vs_hc) = "active_vs_hc"
#Add trait to MEs
ME_AvsHC = orderMEs(cbind(MEs, active_vs_hc))
# Plot the relationships among the eigengenes and the trait
# Plot the dendrogram
pdf(file = "active_vs_hc trait/active_vs_hc - Eigengene dendrogram.pdf", 
    height = 4, width = 8)
par(cex = 1.0)
plotEigengeneNetworks(ME_AvsHC, "Active vs HC - Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf(file = "active_vs_hc trait/active_vs_hc - Eigengene adjacency heatmap.pdf", 
    height = 6, width = 7)
par(cex = 1.0)
plotEigengeneNetworks(ME_AvsHC, "Eigengene adjacency heatmap", marHeatmap = c(7,7,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

#Visualising network of eigengenes comparing active_vs_remission
active_vs_remiss = as.data.frame(cleantraits$active_vs_remission)
names(active_vs_remiss) = "active_vs_remission"
#Add trait to MEs
ME_AvsRM = orderMEs(cbind(MEs, active_vs_remiss))
# Plot the relationships among the eigengenes and the trait
# Plot the dendrogram
pdf(file = "active_vs_remission trait/Active_vs_remission - Eigengene dendrogram.pdf", 
    width = 8, height = 4)
par(cex = 1.0)
plotEigengeneNetworks(ME_AvsRM, "Active vs Remission - Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf(file = "active_vs_remission trait/Active_vs_remission - Eigengene adjacency heatmap.pdf", 
    width = 7, height = 5)
par(cex = 1.0)
plotEigengeneNetworks(ME_AvsRM, "Active vs Remission - Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

#Visualising network of eigengenes comparing remission_vs_hc
remission_vs_hc = as.data.frame(cleantraits$remission_vs_hc)
names(remission_vs_hc) = "remission_vs_hc"
#Add trait to MEs
ME_RMvsHC = orderMEs(cbind(MEs, remission_vs_hc))
# Plot the relationships among the eigengenes and the trait
# Plot the dendrogram
pdf(file = "remission_vs_hc trait/Remission_vs_hc - Eigengene dendrogram.pdf", 
    width = 8, height = 4)
par(cex = 1.0)
plotEigengeneNetworks(ME_RMvsHC, "Remission vs HC - Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf(file = "remission_vs_hc trait/Remission_vs_hc - Eigengene adjacency heatmap.pdf", 
    width = 7, height = 5)
par(cex = 1.0)
plotEigengeneNetworks(ME_RMvsHC, "Remission vs HC - Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


##Comparing eigengene

#Sort phenodata by timepoint
phenodata <- read.csv("Phenodata.csv")
dim(phenodata)
names(phenodata)
phenodata <- as.data.frame(phenodata)
timepointOrder <- c("untimed", "zero", "three", "twelve")
phenodata_ordered <- phenodata %>% 
        arrange(factor(timepoint_short, levels = timepointOrder))
sortedFilename <- phenodata_ordered$filename

#Sort matrix by timepoint
cleanMatrix.sorted <- cleanMatrix[match(sortedFilename, rownames(cleanMatrix)), ]

#Sort MEs by timepoint
MEs.timepoint <- MEs[match(sortedFilename, rownames(MEs)), ]

#Run Plots for magenta module
which.module = "magenta"
ME=MEs.timepoint[, paste("ME",which.module, sep="")]
pdf("Magenta Module/Heatmap and Eigengene expression.pdf",
    width = 16, height = 11)
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanMatrix.sorted[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

pdf("Magenta Module/Eigengene expression timepoint split.pdf",
    width = 16, height = 8)
ggplot(data.frame(ME),aes(seq_along(ME),ME)) +
        geom_bar(stat="identity") +
        geom_vline(xintercept = c(64.5, 132.5, 174.5), col = "red") +
        labs(title = "Purple", x = "Array Sample", y = "Eigengene Expression") +
        theme(plot.title = element_text(color="black", size=14, face="bold"),
              plot.title.position = "plot",
              axis.title.x = element_text(color="black", size=10),
              axis.title.y = element_text(color="black", size=10))
dev.off()


##Hierarchical clustering of samples and heatmap of genes in magenta module
#Create a gene matrix with only genes in magenta module
#geneInfo_active_vs_hc <- read.csv("active_vs_hc trait/active_vs_hc - geneInfo.csv")
#magentaSymbols <- geneInfo_active_vs_hc %>% 
#    subset(moduleColor == "magenta")
#magentaSymbols <- select(magentaSymbols, geneSymbol)
#magentaSymbols <- deframe(magentaSymbols)

#probeToSymbols <- finalAnnots
#rownames(probeToSymbols) <- probeToSymbols$SYMBOL
#probeToSymbols <- probeToSymbols[magentaSymbols, ]
#magentaProbes <- probeToSymbols$PROBEID

#magentaMatrix <- vascFinalMatrix[magentaProbes, ]
#magentaMatrix <- as.data.frame(magentaMatrix)
#magentaMatrix <- t(magentaMatrix)

#Create a gene matrix with only AAV samples
#timeZero <- phenodata[phenodata$timepoint_short == "zero", ]
#timeZeroFile <- timeZero$filename

#magentaMatrix <- magentaMatrix[timeZeroFile, ]
#magentaCluster <- hclust(dist(magentaMatrix), method = "complete")
#plot(magentaCluster, labels = FALSE)

#library(ComplexHeatmap)
#Heatmap(magentaMatrix, show_row_names = FALSE, show_column_names = FALSE)

#Did not show a different between sample clusters


##Exporting purple network to VisANT
load(file = "Vascmicroarray-final.RData")
load(file = "VascNetworkConstruction.RData")
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(cleanMatrix, power = 9);
# Read in the annotation file
finalAnnots
# Select module
module = "magenta";
# Select module probes
probes = colnames(cleanMatrix)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

save(modTOM, TOM, cleanMatrix, file = "Network Analysis.RData")
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0.03,
                            probeToGene = data.frame(finalAnnots$PROBEID, finalAnnots$SYMBOL) )

vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, "-04thresh.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0.04,
                            probeToGene = data.frame(finalAnnots$PROBEID, finalAnnots$SYMBOL) )

#Restrict genes to top 30 hub genes
nTop = 30;
IMConn = softConnectivity(cleanMatrix[, modProbes],
                          power = 9, type = "signed");
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0.03,
                            probeToGene = data.frame(finalAnnots$PROBEID, finalAnnots$SYMBOL) )


##Cytoscape
# Export the network into edge and node list files Cytoscape can read
modGenes = finalAnnots$SYMBOL[match(modProbes, finalAnnots$PROBEID)]
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.04,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

#To onclude all genes, threshold set to 0
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), "allgenes.txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), "allgenes.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.01,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);


##PINK
module = "pink";
# Select module probes
probes = colnames(cleanMatrix)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
modGenes = finalAnnots$SYMBOL[match(modProbes, finalAnnots$PROBEID)]
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("Pink Module/Network/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("Pink Module/Network/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.08,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);
   


#Finding hub genes
chooseTopHubInEachModule(
    cleanMatrix, 
    moduleColors, 
    omitColors = "grey", 
    power = 9, 
    type = "signed")

#library(igraph)
#load("VascNetworkConstruction.RData")
#wgcna2igraph(net, cleanMatrix, top.n.edges = NA, modules2plot = c("magenta","pink"),
#             colors2plot = NULL, kME.threshold = 0.75, adjacency.threshold = 0.1,
#             adj.power = 6, verbose = T, min.edge = 2, node.size = 0,
#             frame.color = NA, node.color = NA, edge.alpha = 0.5, edge.width = 1,
#             returnNet = TRUE)





##TOM plot to visual network
#take out grey genes to save time
#power=9
#color1=mergedColors[net$blockGenes[[1]]]
#restGenes= (color1 != "grey")
#diss1=1-TOMsimilarityFromExpr(cleanMatrix[, restGenes], power = 9 )
#hier1=hclust(as.dist(diss1), method="average" )
#diag(diss1) = NA;
#sizeGrWindow(7,7)
#TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
#        main = "TOM heatmap plot, module genes" )
##Possibly not possible due to number of genes.


##Choosing hub genes in magenta module
load(file = "Vascmicroarray-final.RData")
load(file = "VascNetworkConstruction.RData")
load("Network Analysis.RData")
geneInfo_active_vs_hc <- read.csv("active_vs_hc trait/active_vs_hc - geneInfo.csv")
geneInfoOrdered <- geneInfo_active_vs_hc[order(geneInfo_active_vs_hc$probeID),]
intCon <- intramodularConnectivity(TOM, moduleColors, scaleByMax = FALSE)
IntModCon <- intCon$kWithin
TotalCon <- intCon$kTotal
geneInfo <- cbind(geneInfoOrdered, IntModCon, TotalCon)

magentaInfo <- geneInfo[grep("magenta", geneInfo$moduleColor),]
magentaInfo <- magentaInfo %>% 
    select(probeID, geneSymbol, geneName, GS.active_vs_hc, p.GS.active_vs_hc, MM.magenta, p.MM.magenta, IntModCon, TotalCon)

AbsMM.magenta <- abs(magentaInfo$MM.magenta)
AbsGS <- abs(magentaInfo$GS.active_vs_hc)
magentaInfo <- cbind(magentaInfo, AbsMM.magenta, AbsGS)

magentaHub <- magentaInfo %>%
    subset(AbsGS > 0.5 & AbsMM.magenta > 0.8)

top10MM <- top_n(magentaInfo, 10, magentaInfo$AbsMM.magenta)
top10Con <- top_n(magentaInfo, 10, magentaInfo$IntModCon)

hub <- intersect(top10MM, top10Con)
orderedHub <- hub[order(hub$IntModCon),]
write.csv(orderedHub, "MagentaHubgenes.csv")

orderedHub <- read.csv("MagentaHubgenes.csv")
magenta_nodes <- read.delim("CytoscapeInput-nodes-magenta.txt", header = TRUE, sep = "\t")
in_hub <- (magenta_nodes$altName %in% orderedHub$geneSymbol)
in_node <- magentaInfo$geneSymbol %in% magenta_nodes$altName
in_node <- magentaInfo[in_node,]

magenta_pos <- read.csv("Magenta Module/Magenta Positive Hallmark.csv")
magenta_pos <- select(magenta_pos, Term, Genes)
magenta_neg <- read.csv("Magenta Module/Magenta Negative Hallmark.csv")
magenta_neg <- select(magenta_neg, Term, Genes)
magenta_combine <- rbind(magenta_pos, magenta_neg)
magenta_combine <- magenta_combine %>% 
    mutate(Genes = strsplit(as.character(Genes), ";")) %>% 
    unnest(Genes)
magenta_wide <- magenta_combine %>% 
    group_by(Genes) %>%
    mutate(Term = paste(Term, collapse = ", ")) %>%
    ungroup()
magenta_wide <- as.data.frame(unique(magenta_wide))
rownames(magenta_wide) <- magenta_wide$Genes
magenta_pathways <- magenta_wide[in_node$geneSymbol,]
magenta_final_nodes <- cbind(in_node, in_hub, magenta_pathways)

write.table(magenta_final_nodes, file = "magenta_final_nodes.txt", sep = "\t",
            row.names = FALSE)
write.xlsx(magenta_final_nodes, file = "magenta_final_nodes.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)

#Create pie charts for each gene for images in cytoscape
value <- rep(1, nrow(magenta_combine))
magenta_combine <- cbind(magenta_combine, value)
magenta_combine <- magenta_combine %>%
    group_by(Genes, Term) %>%
    summarise(count = n()) %>%
    group_by(Genes) %>%
    mutate(per=count/sum(count)) %>% 
    ungroup()

nodes <- magenta_final_nodes$geneSymbol
magenta_nodes <- magenta_combine[magenta_combine$Genes %in% nodes,]
magenta_nodes$Term <- str_replace_all(magenta_nodes$Term, 'Pp', 'P')

colourCount = length(unique(magenta_nodes$Term))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

ggplot(magenta_nodes, aes(x="", y=per, fill=Term)) +
    geom_bar(width = 1, stat = "identity") +
    facet_wrap(~Genes) +
    coord_polar("y", start=0) +
    theme_void() +
    theme(axis.text.x=element_blank(), legend.position = "bottom") +
    scale_fill_manual(values = getPalette(colourCount))
ggsave("Magenta Module/Network/Magenta node images.png", height = 15, width = 15)

##Pink hubgenes
pinkInfo <- geneInfo[grep("pink", geneInfo$moduleColor),]
pinkInfo <- pinkInfo %>% 
    select(probeID, geneSymbol, geneName, GS.active_vs_hc, p.GS.active_vs_hc, MM.pink, p.MM.pink, IntModCon, TotalCon)

AbsMM.pink <- abs(pinkInfo$MM.pink)
AbsGS <- abs(pinkInfo$GS.active_vs_hc)
pinkInfo <- cbind(pinkInfo, AbsMM.pink, AbsGS)

pinkHub <- pinkInfo %>%
    subset(AbsGS > 0.5 & AbsMM.pink > 0.8)

top10MM <- top_n(pinkInfo, 10, pinkInfo$AbsMM.pink)
top10Con <- top_n(pinkInfo, 10, pinkInfo$IntModCon)

hub <- intersect(top10MM, top10Con)
orderedHub <- hub[order(hub$IntModCon),]
write.csv(orderedHub, "PinkHubgenes.csv")

orderedHub <- read.csv("PinkHubgenes.csv")
pink_nodes <- read.delim("Pink Module/Network/CytoscapeInput-nodes-pink.txt", header = TRUE, sep = "\t")
in_hub <- (pink_nodes$altName %in% orderedHub$geneSymbol)
in_node <- pinkInfo$geneSymbol %in% pink_nodes$altName
in_node <- pinkInfo[in_node,]

pink_pos <- read.csv("Pink Module/Pink Positive Hallmark.csv")
pink_pos <- select(pink_pos, Term, Genes)
pink_neg <- read.csv("Pink Module/Pink Negative Hallmark.csv")
pink_neg <- select(pink_neg, Term, Genes)
pink_combine <- rbind(pink_pos, pink_neg)
pink_combine <- pink_combine %>% 
    mutate(Genes = strsplit(as.character(Genes), ";")) %>% 
    unnest(Genes)
pink_wide <- pink_combine %>% 
    group_by(Genes) %>%
    mutate(Term = paste(Term, collapse = ", ")) %>%
    ungroup()
pink_wide <- as.data.frame(unique(pink_wide))
rownames(pink_wide) <- pink_wide$Genes
pink_pathways <- pink_wide[in_node$geneSymbol,]
pink_final_nodes <- cbind(in_node, in_hub, pink_pathways)

write.xlsx(pink_final_nodes, file = "Pink Module/Network/pink_final_nodes.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)

#Create pie charts for each gene for images in cytoscape
value <- rep(1, nrow(pink_combine))
pink_combine <- cbind(pink_combine, value)
pink_combine <- pink_combine %>%
    group_by(Genes, Term) %>%
    summarise(count = n()) %>%
    group_by(Genes) %>%
    mutate(per=count/sum(count)) %>% 
    ungroup()

nodes <- pink_final_nodes$geneSymbol
pink_nodes <- pink_combine[pink_combine$Genes %in% nodes,]
pink_nodes$Term <- str_replace_all(pink_nodes$Term, 'Pp', 'P')

colourCount = length(unique(pink_nodes$Term))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

ggplot(pink_nodes, aes(x="", y=per, fill=Term)) +
    geom_bar(width = 1, stat = "identity") +
    facet_wrap(~Genes) +
    coord_polar("y", start=0) +
    theme_void() +
    theme(axis.text.x=element_blank(), legend.position = "bottom") +
    scale_fill_manual(values = getPalette(colourCount))
ggsave("Pink Module/Network/Pink node images.png", height = 15, width = 15)
