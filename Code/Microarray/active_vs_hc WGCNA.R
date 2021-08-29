##WGCNA
#active_vs_hc trait module analysis

library(WGCNA)
library(fgsea)
setwd("~/VascMicroarray/Usable Data")

options(stringsAsFactors = FALSE);

wcgnaAnalysis = load(file = "VascNetworkConstruction.RData")
wcgnaAnalysis = load(file = "VascWGCNA_DataInput.RData")

# Define numbers of genes and samples
nGenes = ncol(cleanMatrix);
nSamples = nrow(cleanMatrix);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(cleanMatrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, cleantraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#active_vs_hc is  significant in magenta, pink modules

active_vs_hc = as.data.frame(cleantraits$active_vs_hc)
names(active_vs_hc) = "active_vs_hc"

#Names(colours) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(cleanMatrix, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance_active_vs_hc = as.data.frame(cor(cleanMatrix, active_vs_hc, use = "p"));
GSPvalue_active_vs_hc = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_active_vs_hc), 
                                                              nSamples));
names(geneTraitSignificance_active_vs_hc) = paste("GS.", names(active_vs_hc), sep="");
names(GSPvalue_active_vs_hc) = paste("p.GS.", names(active_vs_hc), sep="");

#Plot Module Significance in active_vs_hc
pdf(file = "active_vs_hc trait/active_vs_hc - GS across modules.pdf",
    width = 8, height = 4)
plotModuleSignificance(t(abs(geneTraitSignificance_active_vs_hc)),
                       moduleColors,
                       boxplot = FALSE,
                       main = "Gene significance across modules in Active vs HC trait",
                       ylab = "Gene significance",
                       cex.names = 0.8,
                       las = 2)
dev.off()

#Intramodular analysis: identifying genes with high GS and MM
mm.gs <- cbind(geneModuleMembership, geneTraitSignificance_active_vs_hc, moduleColors)
magenta.mm.gs <- mm.gs[grep("magenta", mm.gs$moduleColors),]
pink.mm.gs <- mm.gs[grep("pink", mm.gs$moduleColors),]

pdf(file = "active_vs_hc trait/active_vs_hc - magenta - MM-GS.pdf", 
    width = 8, height = 4)
ggplot(magenta.mm.gs, aes(x = abs(MMmagenta), y = abs(GS.active_vs_hc))) +
  geom_point(color = "magenta") +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "Module Membership in magenta module", y = "Gene significance for Active vs HC") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 0.8) +
  xlim(0, 1) +
  ylim(0, 1)
dev.off()

pdf(file = "active_vs_hc trait/active_vs_hc - pink - MM-GS.pdf", 
    width = 8, height = 4)
ggplot(pink.mm.gs, aes(x = abs(MMpink), y = abs(GS.active_vs_hc))) +
  geom_point(color = "pink") +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "Module Membership in pink module", y = "Gene significance for Active vs HC") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 0.8) +
  xlim(0, 1) +
  ylim(0, 1)
dev.off()



pdf(file = "active_vs_hc trait/active_vs_hc - pink - MM-GS.pdf", 
    width = 8, height = 4)
module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance_active_vs_hc[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Active vs HC",
                   main = paste("Module membership vs. Gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                   col = module, pch = 19)
dev.off()


#Load in gene annotations
probes = colnames(cleanMatrix)
probes2annot = match(probes, finalAnnots$PROBEID)
#cleanExpMatrix = wgcna_expmatrix[, probes2annot]
#probes = names(cleanExpMatrix)
sum(is.na(probes2annot))


##active_vs_hc - Create a table with genes associated with modules
# Create the starting data frame
geneInfo0_active_vs_hc = data.frame(probeID = probes,
                                           geneSymbol = finalAnnots$SYMBOL[probes2annot],
                                           geneName = finalAnnots$GENENAME[probes2annot],
                                           moduleColor = moduleColors,
                                           geneTraitSignificance_active_vs_hc,
                                           GSPvalue_active_vs_hc)

# Order modules by their significance for active_vs_hc
modOrder = order(-abs(cor(MEs, active_vs_hc, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0_active_vs_hc)
  geneInfo0_active_vs_hc = data.frame(geneInfo0_active_vs_hc, geneModuleMembership[, modOrder[mod]],
                                             MMPvalue[, modOrder[mod]]);
  names(geneInfo0_active_vs_hc) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                           paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_active_vs_hc = order(geneInfo0_active_vs_hc$moduleColor, 
                                      -abs(geneInfo0_active_vs_hc$GS.active_vs_hc));
geneInfo_active_vs_hc = geneInfo0_active_vs_hc[geneOrder_active_vs_hc, ]

#save
write.csv(geneInfo_active_vs_hc, 
          file = "active_vs_hc trait/active_vs_hc - geneInfo.csv")

#FGSEA
library(fgsea)
library(tidyverse)

#Load in DEG for active_vs_hc
#active_vs_hc_deg <- read.csv("active_vs_hc trait/active_vs_hc_deg.csv")
#active_vs_hc_deg <- active_vs_hc_deg %>% 
#  select(SYMBOL, logFC)
#rownames(active_vs_hc_deg) <- active_vs_hc_deg$SYMBOL

#Create custom gene set 
setwd("Gene sets/")
geneList = list.files(full.names = TRUE, recursive = T, patt ="*.txt")
geneListFinal = lapply(geneList, read.delim)
geneListFinal = unlist(geneListFinal, recursive = FALSE)
setwd("~/VascMicroarray/Usable Data")

pathways.hallmark <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

#FGSEA of genes in magenta module
#Get symbols and module membership for magenta module
geneInfo_active_vs_hc <- read.csv("active_vs_hc trait/active_vs_hc - geneInfo.csv")
magentaMM <- geneInfo_active_vs_hc %>% 
  subset(moduleColor == "magenta")
magentaMM <- select(magentaMM, geneSymbol, MM.magenta)
magentaMM <- deframe(magentaMM)

#Run FGSEA in magenta module on custom gene set
magentaFgseaMM_custom <- fgsea(pathways = geneListFinal, stats = magentaMM, nperm = 1000)
magentaFgseaMM_custom <- magentaFgseaMM_custom[order(pval), ]
magentaFgseaMM_custom <- as.data.frame(magentaFgseaAvsHC_custom)
magentaFgseaMM_custom <- apply(magentaFgseaMM_custom,2,as.character)
write.csv(magentaFgseaMM_custom, file = "Magenta Module/Magenta module membership FGSEA - custom.csv")

topPathwaysUp_custom <- magentaFgseaMM_custom[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- magentaFgseaMM_custom[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Magenta Module/ Magenta module membership GSEA table - custom.pdf", 
    height = 5, width = 10)
plotGseaTable(geneListFinal[topPathways_custom], magentaMM, magentaFgseaMM_custom, 
              gseaParam = 0.5)
dev.off()

plotEnrichment(geneListFinal[["GO_MACROPHAGE_ACTIVATION"]],
               magentaMM) + 
  labs(title = "GO_MACROPHAGE_ACTIVATION")

#Run FGSEA in magenta module on hallmark gene set
magentaFgseaMM_hallmark <- fgsea(pathways = pathways.hallmark, stats = magentaMM, nperm = 1000)
magentaFgseaMM_hallmark <- magentaFgseaMM_hallmark[order(pval), ]
magentaFgseaMM_hallmark <- as.data.frame(magentaFgseaMM_hallmark)
magentaFgseaMM_hallmark <- apply(magentaFgseaMM_hallmark,2,as.character)
write.csv(magentaFgseaMM_hallmark, file = "Magenta Module/Magenta module membership FGSEA - hallmark.csv")

topPathwaysUp_hallmark <- magentaFgseaMM_hallmark[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_hallmark <- magentaFgseaMM_hallmark[ES < 0][head(order(pval), n=10), pathway]
topPathways_hallmark <- c(topPathwaysUp_hallmark, rev(topPathwaysDown_hallmark))
pdf(file = "Magenta Module/ Magenta module membership GSEA table - hallmark.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_hallmark], magentaMM, magentaFgseaMM_hallmark, 
              gseaParam = 0.5)
dev.off()

pdf(file = "Magenta Module/HALLMARK_COMPLEMENT FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_COMPLEMENT"]],
               magentaMM) + 
  labs(title = "HALLMARK_COMPLEMENT")
dev.off()

pdf(file = "Magenta Module/HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]],
               magentaMM) + 
  labs(title = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
dev.off()

pdf(file = "Magenta Module/HALLMARK_OXIDATIVE_PHOSPHORYLATION FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
               magentaMM) + 
  labs(title = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
dev.off()

pdf(file = "Magenta Module/HALLMARK_INFLAMMATORY_RESPONSE FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
               magentaMM) + 
  labs(title = "HALLMARK_INFLAMMATORY_RESPONSE")
dev.off()

pdf(file = "Magenta Module/HALLMARK_IL2_STAT5_SIGNALING FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]],
               magentaMM) + 
  labs(title = "HALLMARK_IL2_STAT5_SIGNALING")
dev.off()

pdf(file = "Magenta Module/HALLMARK_IL6_JAK_STAT3_SIGNALING FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]],
               magentaMM) + 
  labs(title = "HALLMARK_IL6_JAK_STAT3_SIGNALING")
dev.off()

#FGSEA of genes in pink module
#Get symbols and module membership for pink module
geneInfo_active_vs_hc <- read.csv("active_vs_hc trait/active_vs_hc - geneInfo.csv")
pinkMM <- geneInfo_active_vs_hc %>% 
  subset(moduleColor == "pink")
pinkMM <- select(pinkMM, geneSymbol, MM.pink)
pinkMM <- deframe(pinkMM)

#Run FGSEA in pink module on custom gene set
pinkFgseaMM_custom <- fgsea(pathways = geneListFinal, stats = pinkMM, nperm = 1000)
pinkFgseaMM_custom <- pinkFgseaMM_custom[order(pval), ]
pinkFgseaMM_custom <- as.data.frame(pinkFgseaMM_custom)
pinkFgseaMM_custom <- apply(pinkFgseaMM_custom,2,as.character)
write.csv(pinkFgseaMM_custom, file = "Pink module/Pink module membership FGSEA - custom.csv")

topPathwaysUp_custom <- pinkFgseaMM_custom[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- pinkFgseaMM_custom[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "Pink module/Pink module membership GSEA table - custom.pdf", 
    height = 5, width = 10)
plotGseaTable(geneListFinal[topPathways_custom], pinkMM, pinkFgseaMM_custom, 
              gseaParam = 0.5)
dev.off()

#Get symbols and logFC for pink module
pinkGenes <- read.csv("Pink module/pinkGenes.csv")
pinkGenes <- pinkGenes$x
deg_active_vs_hc <- read.csv("active_vs_hc trait/active_vs_hc_deg.csv")
rownames(deg_active_vs_hc) <- deg_active_vs_hc$SYMBOL
pink_logfc <- deg_active_vs_hc[pinkGenes,]
pink_logfc <- pink_logfc %>% 
  select(SYMBOL, logFC)
pink_logfc <- deframe(pink_logfc)

#Run FGSEA in pink module on custom gene set
pinkFgseaAvsHC_custom <- fgsea(pathways = geneListFinal, stats = pink_logfc, nperm = 1000)
pinkFgseaAvsHC_custom <- pinkFgseaAvsHC_custom[order(pval), ]
pinkFgseaAvsHC_custom <- as.data.frame(pinkFgseaAvsHC_custom)
pinkFgseaAvsHC_custom <- apply(pinkFgseaAvsHC_custom,2,as.character)
write.csv(pinkFgseaAvsHC_custom, file = "active_vs_hc trait/active_vs_hc - pink module FGSEA - custom.csv")

topPathwaysUp_custom <- pinkFgseaAvsHC_custom[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_custom <- pinkFgseaAvsHC_custom[ES < 0][head(order(pval), n=10), pathway]
topPathways_custom <- c(topPathwaysUp_custom, rev(topPathwaysDown_custom))
pdf(file = "active_vs_hc trait/active_vs_hc - pink module GSEA table - custom.pdf", 
    height = 5, width = 10)
plotGseaTable(geneListFinal[topPathways_custom], pink_logfc, pinkFgseaAvsHC_custom, 
              gseaParam = 0.5)
dev.off()


#Run FGSEA in pink module on hallmark gene set
pinkFgseaMM_hallmark <- fgsea(pathways = pathways.hallmark, stats = pinkMM, nperm = 1000)
pinkFgseaMM_hallmark <- pinkFgseaMM_hallmark[order(pval), ]
pinkFgseaMM_hallmark <- as.data.frame(pinkFgseaMM_hallmark)
pinkFgseaMM_hallmark <- apply(pinkFgseaMM_hallmark,2,as.character)
write.csv(pinkFgseaMM_hallmark, file = "Pink module/ Pink module membership FGSEA - hallmark.csv")

topPathwaysUp_hallmark <- pinkFgseaMM_hallmark[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_hallmark <- pinkFgseaMM_hallmark[ES < 0][head(order(pval), n=10), pathway]
topPathways_hallmark <- c(topPathwaysUp_hallmark, rev(topPathwaysDown_hallmark))
pdf(file = "Pink module/ Pink module membership GSEA table - hallmark.pdf", 
    height = 5, width = 10)
plotGseaTable(pathways.hallmark[topPathways_hallmark], pinkMM, pinkFgseaMM_hallmark, 
              gseaParam = 0.5)
dev.off()

pdf(file = "Pink module/HALLMARK_IL6_JAK_STAT3_SIGNALING FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]],
               pinkMM) + 
  labs(title = "HALLMARK_IL6_JAK_STAT3_SIGNALING")
dev.off()

pdf(file = "Pink module/HALLMARK_IL2_STAT5_SIGNALING FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]],
               pinkMM) + 
  labs(title = "HALLMARK_IL2_STAT5_SIGNALING")
dev.off()

pdf(file = "Pink module/HALLMARK_INFLAMMATORY_RESPONSE FGSEA MM.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
               pinkMM) + 
  labs(title = "HALLMARK_INFLAMMATORY_RESPONSE")
dev.off()

pdf(file = "Pink module/HALLMARK_P53_PATHWAY.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_P53_PATHWAY"]],
               pinkMM) + 
  labs(title = "HALLMARK_P53_PATHWAY")
dev.off()

pdf(file = "Pink module/HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf", 
    height = 3, width = 6)
plotEnrichment(pathways.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               pinkMM) + 
  labs(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
dev.off()

