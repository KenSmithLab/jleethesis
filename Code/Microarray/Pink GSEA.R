##Pink GSEA

setwd("~/VascMicroarray/Usable Data")
library(tidyverse)

##Split upregated genes and downregulated genes between active and hc
#Get genes from pink module
pinkGenes <- read.csv("active_vs_hc trait/active_vs_hc - geneInfo.csv")
pinkGenes <- pinkGenes %>% 
  subset(moduleColor == "pink")
pinkGenes <- dplyr::select(pinkGenes, geneSymbol)
pinkGenes <- deframe(pinkGenes)
write.csv(pinkGenes, file = "Pink Module/pinkGenes.csv")

#Get logFC
active_vs_hc_deg <- read.csv("active_vs_hc trait/active_vs_hc_deg.csv")
rownames(active_vs_hc_deg) <- active_vs_hc_deg$SYMBOL
pink_logfc <- active_vs_hc_deg[pinkGenes, ]
pink_logfc <- pink_logfc %>% 
  dplyr::select(SYMBOL, logFC)

pink_pos <- pink_logfc %>% 
  subset(logFC > 0)
pink_pos <- as.vector(pink_pos$SYMBOL)
pink_neg <- pink_logfc %>% 
  subset(logFC < 0)
pink_neg <- as.vector(pink_neg$SYMBOL)


##EnrichR for gene set enrichment in pink module
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Biological_Process_2018", "KEGG_2019_Human", "Reactome_2016", "MSigDB_Hallmark_2020")

#pink Positive genes
if (websiteLive) {
  enriched_pinkpos <- enrichr(pink_pos, dbs)
}

go_bp_pinkpos <- if (websiteLive) enriched_pinkpos[["GO_Biological_Process_2018"]]
write.csv(go_bp_pinkpos, "Pink module/Pink Positive GO_Biological_Process.csv")
pdf(file = "Pink module/Pink Positive GO_Biological_Process.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkpos[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

kegg_pinkpos <- if (websiteLive) enriched_pinkpos[["KEGG_2019_Human"]]
write.csv(kegg_pinkpos, "Pink module/Pink Positive KEGG.csv")
pdf(file = "Pink module/Pink Positive KEGG.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkpos[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

reactome_pinkpos <- if (websiteLive) enriched_pinkpos[["Reactome_2016"]]
write.csv(reactome_pinkpos, "Pink module/Pink Positive Reactome.csv")
pdf(file = "Pink module/Pink Positive Reactome.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkpos[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

hallmark_pinkpos <- if (websiteLive) enriched_pinkpos[["MSigDB_Hallmark_2020"]]
write.csv(hallmark_pinkpos, "Pink module/Pink Positive Hallmark.csv")
pdf(file = "Pink module/Pink Positive Hallmark.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkpos[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#pink Negative genes
if (websiteLive) {
  enriched_pinkneg <- enrichr(pink_neg, dbs)
}

go_bp_pinkneg <- if (websiteLive) enriched_pinkneg[["GO_Biological_Process_2018"]]
write.csv(go_bp_pinkneg, "Pink module/Pink Negative GO_Biological_Process.csv")
pdf(file = "Pink module/Pink Negative GO_Biological_Process.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkneg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

kegg_pinkneg <- if (websiteLive) enriched_pinkneg[["KEGG_2019_Human"]]
write.csv(kegg_pinkneg, "Pink module/Pink Negative KEGG.csv")
pdf(file = "Pink module/Pink Negative KEGG.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkneg[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

reactome_pinkneg <- if (websiteLive) enriched_pinkneg[["Reactome_2016"]]
write.csv(reactome_pinkneg, "Pink module/Pink Negative Reactome.csv")
pdf(file = "Pink module/Pink Negative Reactome.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkneg[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

hallmark_pinkneg <- if (websiteLive) enriched_pinkneg[["MSigDB_Hallmark_2020"]]
write.csv(hallmark_pinkneg, "Pink Module/Pink Negative Hallmark.csv")
pdf(file = "Pink module/Pink Negative Hallmark.pdf", width = 8, height = 5)
if (websiteLive) plotEnrich(enriched_pinkneg[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()


##Genes from EnrichR Pathways
load(file = "Vascmicroarray-final.RData")
phenodata <- read.csv("Phenodata.csv")
phenodata <- as.data.frame(phenodata)
timepoint <- as.character(phenodata$timepoint_short)
diagnosis <- as.character(phenodata$diagnosis_short)
study_id <- as.character(phenodata$study_id)
anca <- as.character(phenodata$ANCA_short)

active_vs_hc_deg <- read.csv("active_vs_hc trait/active_vs_hc_deg.csv")
rownames(active_vs_hc_deg) <- active_vs_hc_deg$SYMBOL

genes <- c("PRDX2", "PRDX6", "SOD1", "SOD2", "NDUFB4", 
           "ALAS1", "MPC1", "UQCR11", "COX6A1", "ATP6V1E1", 
           "IDH3A", "ITGAM", "C1R", "GRB2", "LTA4H", 
           "IL15", "IL18", "PROK2", "KCNJ2", "C5AR1")
genesAnnots <- active_vs_hc_deg[genes,]
genesAnnots <- genesAnnots[, c("PROBEID", "GENENAME", "SYMBOL")]
genesProbes <- as.character(genesAnnots$PROBEID)

genesExpr <- as.data.frame(vascFinalMatrix[genesProbes,])
rownames(genesExpr) <- genesAnnots$SYMBOL
genesExpr <- as.data.frame(t(genesExpr))
genesExpr <- cbind(genesExpr, timepoint, diagnosis, study_id, anca)

genesLong <- genesExpr %>% 
  pivot_longer(cols = PRDX2:C5AR1,
               names_to = "gene",
               values_to = "value")
genesLong$timepoint <- factor(genesLong$timepoint,levels = c("untimed", "zero", "three", "twelve"))
genesLong$gene <- factor(genesLong$gene,levels = c("PRDX2", "PRDX6", "SOD1", "SOD2", "NDUFB4", 
                                                   "ALAS1", "MPC1", "UQCR11", "COX6A1", "ATP6V1E1", 
                                                   "IDH3A", "ITGAM", "C1R", "GRB2", "LTA4H", 
                                                   "IL15", "IL18", "PROK2", "KCNJ2", "C5AR1"))

pdf(file = "pink Module/Expression across timepoints for selected genes.pdf", width = 15, height = 15)
ggplot(genesLong, aes(x = gene, y = value, fill = timepoint)) +
  geom_boxplot() +
  facet_wrap(~gene, scale="free")
dev.off()

genesRos <- genesExpr %>% 
  select(PRDX2, PRDX6, SOD1, SOD2, NDUFB4, timepoint, diagnosis, study_id, anca)
genesRos$timepoint <- factor(genesRos$timepoint,levels = c("untimed", "zero", "three", "twelve"))

genesRosLong <- genesRos %>% 
  pivot_longer(cols = PRDX2:NDUFB4,
               names_to = "gene",
               values_to = "value")
genesRosLong$timepoint <- factor(genesRosLong$timepoint,
                                 levels = c("untimed", "zero", "three", "twelve"))

pdf(file = "pink Module/ROS genes expression.pdf", width = 8, height = 6)
ggplot(genesRosLong, aes(x = gene, y = value, fill = timepoint)) +
  geom_boxplot() +
  facet_wrap(~gene, scale="free") +
  labs(x = "Gene", y = "Expression", title = "Reactive Oxidative Species Pathway") +
  theme_classic() +
  scale_fill_manual(values = c("skyblue1", "brown1", "darkorange1", "yellow2"),
                    name = "Timepoint (months)",
                    labels = c("Healthy", "AAV Zero", "AAV Three", "AAV Twelve")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#NDUFB4 gene split by PR3 and MPO
split_genesRos <- genesRos[diagnosis != "Healthy",]
split_genesRos <- subset(split_genesRos, anca != "double")
split_genesRos <- subset(split_genesRos, anca != "neg")
split_genesRos$timepoint <- factor(split_genesRos$timepoint,levels = c("zero", "three", "twelve"))
split_genesRos$anca <- factor(split_genesRos$anca, levels = c("PR3", "MPO"))

#Stats between pink and ANCA intersecting with Time
split_ndufb4_aov <- aov(NDUFB4 ~ anca * timepoint, data = split_genesRos)
# Summary of the analysis
summary(split_ndufb4_aov)
#Determine which groups are significant
split_ndufb4_bonf <- split_genesRos %>% 
  group_by(timepoint) %>%
  pairwise_t_test(
    NDUFB4 ~ anca,
    p.adjust.method = "bonferroni"
  ) %>% 
  ungroup()

ndufb4_aov <- aov(NDUFB4 ~ timepoint, data = genesRos)
summary(ndufb4_aov)

ndufb4_bonf <- genesRos %>% 
  pairwise_t_test(NDUFB4 ~ timepoint,
                  p.adjust.method = "bonferroni")
ndufb4_bonf

ggplot(genesRos, aes(x = timepoint, y = NDUFB4, fill = timepoint)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)

#Oxidative Phosphorylation
genesOxPhos <- genesExpr %>% 
  select(ALAS1, MPC1, UQCR11, COX6A1, ATP6V1E1, IDH3A, timepoint, diagnosis, study_id, anca)
genesOxPhosLong <- genesOxPhos %>% 
  pivot_longer(cols = ALAS1:IDH3A,
               names_to = "gene",
               values_to = "value")

pdf(file = "pink Module/OxPhos genes expression.pdf", width = 8, height = 6)
ggplot(genesOxPhosLong, aes(x = gene, y = value, fill = timepoint)) +
  geom_boxplot() +
  facet_wrap(~gene, scale="free") +
  labs(x = "Gene", y = "Expression", title = "Oxidative Phosphorylation Pathway") +
  theme_classic() +
  scale_fill_manual(values = c("skyblue1", "brown1", "darkorange1", "yellow2"),
                    name = "Timepoint (months)",
                    labels = c("Healthy", "AAV Zero", "AAV Three", "AAV Twelve")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

##Complement
genesComp <- genesExpr %>% 
  select(ITGAM, C1R, GRB2, LTA4H, timepoint, diagnosis, study_id, anca)
genesCompLong <- genesComp %>% 
  pivot_longer(cols = ITGAM:LTA4H,
               names_to = "gene",
               values_to = "value")

#Stats between pink ME and Diagnosis
comp_aov <- aov(value ~ timepoint * gene, data = genesCompLong)
# Summary of the analysis
summary(comp_aov)

comp_bonf <- genesCompLong %>%
  group_by(gene) %>%
  t_test(value ~ timepoint) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  ungroup()
comp_bonf

comp_bonf <- comp_bonf %>% add_xy_position(x = "timepoint", scales = "free", step.increase = 0.5)

pdf(file = "pink Module/Complement genes expression.pdf", width = 6, height = 5)
ggboxplot(genesCompLong, x = "timepoint", y = "value", fill = "timepoint", facet.by = "gene") +
  facet_wrap(~gene, scales = "free") +
  labs(x = "Gene", y = "Expression", title = "Complement Pathway") +
  theme_classic() +
  scale_fill_manual(values = c("skyblue1", "brown1", "darkorange1", "yellow2"),
                    name = "Timepoint (months)",
                    labels = c("Healthy", "AAV Zero", "AAV Three", "AAV Twelve")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(comp_bonf, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
dev.off()

#Inflammatory Response
genesInflam <- genesExpr %>% 
  select(IL15, IL18, PROK2, KCNJ2, C5AR1, timepoint, diagnosis, study_id, anca)
genesInflamLong <- genesInflam %>% 
  pivot_longer(cols = IL15:C5AR1,
               names_to = "gene",
               values_to = "value")

pdf(file = "pink Module/Inflammatory response genes expression.pdf", width = 8, height = 6)
ggplot(genesInflamLong, aes(x = gene, y = value, fill = timepoint)) +
  geom_boxplot() +
  facet_wrap(~gene, scale="free") +
  labs(x = "Gene", y = "Expression", title = "Inflammatory Response Genes") +
  theme_classic() +
  scale_fill_manual(values = c("skyblue1", "brown1", "darkorange1", "yellow2"),
                    name = "Timepoint (months)",
                    labels = c("Healthy", "AAV Zero", "AAV Three", "AAV Twelve")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()