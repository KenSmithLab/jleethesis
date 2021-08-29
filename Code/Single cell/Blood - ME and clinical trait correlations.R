##Correlating Single cell module eigengene with clinical traits

setwd("~/Jacinta/")

library(WGCNA)
library(ggplot2)
library(tidyr)
library(rstatix)
library(ggpubr)
library(dplyr)

load("~/Jacinta/Blood/Module Preservation/myeloid_all_matrix.RData")

#Get Magenta Eigengenes
magentaGenes <- read.csv("magentaGenes.csv")
magentaGenes <- magentaGenes$x

magenta_myeloid <- myeloid_all_matrix[magentaGenes,]
magenta_myeloid <- t(magenta_myeloid)

magenta <- rep("magenta", times = length(magentaGenes))

magentaME <- moduleEigengenes(magenta_myeloid, magenta)

#Get Pink Eigengenes
pinkGenes <- read.csv("pinkGenes.csv")
pinkGenes <- pinkGenes$x

pink_myeloid <- myeloid_all_matrix[pinkGenes,]
pink_myeloid <- t(pink_myeloid)

pink <- rep("pink", times = length(pinkGenes))

pinkME <- moduleEigengenes(pink_myeloid, pink)

#Load in clinical traits and make a data frame with ME and traits
demo <- read.csv("single_cell_phenodata.csv")

pbmc_demo <- demo[grep("pbmc", demo$filename),]
pbmc_demo <- pbmc_demo[-grep("cv15_pbmc_w12", pbmc_demo$filename),]
MEmagenta <- magentaME$eigengenes$MEmagenta
MEpink <- pinkME$eigengenes$MEpink

demo_ME <- cbind(pbmc_demo, MEpink, MEmagenta)
demo_ME <- as.data.frame(demo_ME)

demo_ME <- transform(demo_ME, crp = as.numeric(crp),
                  bvas_score = as.numeric(bvas_score),
                  gfr = as.numeric(gfr),
                  pr3_titre = as.numeric(pr3_titre))
demo_ME$timepoint <- factor(demo_ME$timepoint, 
                         levels = c("control", "D1", "W12"))

setwd("Blood/ME and Clinical/")
#Stats between MEpink and Time
pink_aov <- aov(MEpink ~ timepoint, data = demo_ME)
# Summary of the analysis
summary(pink_aov)
#Determine which groups are significant
pink_tuk <- demo_ME %>% 
  tukey_hsd(MEpink ~ timepoint)
pink_tuk <- pink_tuk %>% 
  add_xy_position(x = "timepoint")
pink_tuk

pdf(file = "Blood - Pink module over time.pdf", height = 3.5, width = 4)
ggplot(demo_ME, aes(x = timepoint, y = MEpink)) +
  geom_boxplot(outlier.shape = NA, fill = c("skyblue1", "brown1", "darkorange1")) +
  scale_x_discrete(labels = c("HC", "Relapse", "Remission")) +
  geom_jitter(width = 0.1, size = 0.7)+ 
  theme_classic() +
  labs(x = element_blank(), y = "Pink Eigengene") +
  stat_pvalue_manual(pink_tuk, label = "p.adj", tip.length = 0.01)
dev.off()

#Stats between MEmagenta and Time
magenta_aov <- aov(MEmagenta ~ timepoint, data = demo_ME)
# Summary of the analysis
summary(magenta_aov)
#Determine which groups are significant
magenta_tuk <- demo_ME %>% 
  tukey_hsd(MEmagenta ~ timepoint)
magenta_tuk <-magenta_tuk %>% 
  add_xy_position(x = "timepoint")
magenta_tuk

pdf(file = "Blood - Magenta module over time.pdf", height = 3.5, width = 4)
ggplot(demo_ME, aes(x = timepoint, y = MEmagenta)) +
  geom_boxplot(outlier.shape = NA, fill = c("skyblue1", "brown1", "darkorange1")) +
  scale_x_discrete(labels = c("HC", "Relapse", "Remission")) +
  geom_jitter(width = 0.1, size = 0.7)+ 
  theme_classic() +
  labs(x = element_blank(), y = "Magenta Eigengene") +
  ylim(-0.005, 0.0105) +
  stat_pvalue_manual(magenta_tuk, label = "p.adj", tip.length = 0.01)
dev.off()

#CRP and MEpink
crp <- subset(demo_ME, !is.na(crp))

pdf(file = "Blood - Pink module and CRP correlation.pdf", height = 3, width = 5)
ggplot(crp, aes(x = crp, y = MEpink)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint",
                     labels = c("Relapse", "Remission")) +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint",
                    labels = c("Relapse", "Remission")) +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "CRP (mg/dL)", y = "Pink Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 55, label.y = -0.1)
dev.off()

#BVAS and MEpink
bvas <- subset(demo_ME, !is.na(bvas))

pdf(file = "Blood - Pink module and BVAS correlation.pdf", 
    height = 3, width = 5)
ggplot(bvas, aes(x = bvas_score, y = MEpink)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint",
                     labels = c("Relapse", "Remission")) +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint",
                    labels = c("Relapse", "Remission")) +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "BVAS Score", y = "Pink Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 13, label.y = -0.05)
dev.off()

#GFR and MEpink
gfr <- subset(demo_ME, !is.na(gfr))

pdf(file = "Blood - Pink module and GFR correlation.pdf", 
    height = 3, width = 5)
ggplot(gfr, aes(x = gfr, y = MEpink)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint",
                     labels = c("Relapse", "Remission")) +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint",
                    labels = c("Relapse", "Remission")) +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = expression(paste("GFR (mL/min/1.73m" ^ "2", ")")), y = "Pink Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 25, label.y = 0.25)
dev.off()

#pr3_titre and MEpink
pr3 <- subset(demo_ME, !is.na(pr3_titre))

pdf(file = "Blood - Pink module and PR3 titre correlation.pdf", 
    height = 3, width = 5)
ggplot(pr3, aes(x = pr3_titre, y = MEpink)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint") +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint") +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "PR3 Titre (iu/mL)", y = "Pink Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 37, label.y = 0.25)
dev.off()

#CRP and MEmagenta
pdf(file = "Blood - Magenta module and CRP correlation.pdf", height = 3, width = 5)
ggplot(crp, aes(x = crp, y = MEmagenta)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint") +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint") +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "CRP (mg/dL)", y = "Magenta Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 25, label.y = -0.003)
dev.off()

#BVAS and MEmagenta
pdf(file = "Blood - Magenta module and BVAS correlation.pdf", 
    height = 3, width = 5)
ggplot(bvas, aes(x = bvas_score, y = MEmagenta)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint") +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint") +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "BVAS Score", y = "Magenta Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 10, label.y = -0.003)
dev.off()

#GFR and MEmagenta
pdf(file = "Blood - Magenta module and GFR correlation.pdf", 
    height = 3, width = 5)
ggplot(gfr, aes(x = gfr, y = MEmagenta)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint") +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint") +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = expression(paste("GFR (mL/min/1.73m" ^ "2", ")")), y = "Magenta Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 25, label.y = 0.005)
dev.off()

#pr3_titre and MEmagenta
pdf(file = "Blood - Magenta module and PR3 titre correlation.pdf", 
    height = 3, width = 5)
ggplot(pr3, aes(x = pr3_titre, y = MEmagenta)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint") +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint") +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "PR3 Titre (iu/mL)", y = "Magenta Eigengene") +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 37, label.y = 0.005)
dev.off()



##Disease activity in single cell dataset
#Stats between CRP and Time
crp$timepoint <- as.character(crp$timepoint)
crp_t <- crp %>% 
  t_test(crp ~ timepoint)
crp_t

pdf(file = "Blood - CRP levels over time.pdf", height = 3, width = 4)
ggplot(crp, aes(x = timepoint, y = crp)) +
  geom_boxplot(outlier.shape = NA, fill = c("brown1", "darkorange1")) +
  scale_x_discrete(labels = c("Relapse", "Remission")) +
  geom_jitter(width = 0.1, size = 0.7)+ 
  theme_classic() +
  labs(x = element_blank(), y = "CRP (mg/dL)")
dev.off()

#Stats between BVAS and Time
bvas$timepoint <- as.character(bvas$timepoint)
bvas_t <- bvas %>% 
  t_test(bvas_score ~ timepoint) %>% 
  add_significance()
bvas_t

pdf(file = "Blood - BVAS levels over time.pdf", height = 3, width = 4)
ggplot(bvas_magenta, aes(x = timepoint, y = bvas_score)) +
  geom_boxplot(outlier.shape = NA, fill = c("brown1", "darkorange1")) +
  scale_x_discrete(labels = c("Relapse", "Remission")) +
  geom_jitter(width = 0.1, size = 0.7)+ 
  theme_classic() +
  labs(x = element_blank(), y = "BVAS Score") +
  stat_pvalue_manual(bvas_t, label = "p.signif", tip.length = 0.01, y.position = 23)
dev.off()


#Stats between GFR and Time
gfr$timepoint <- as.character(gfr$timepoint)
gfr_t <- gfr %>% 
  t_test(gfr ~ timepoint)
gfr_t

pdf(file = "Blood - GFR levels over time.pdf", height = 3, width = 4)
ggplot(gfr, aes(x = timepoint, y = gfr)) +
  geom_boxplot(outlier.shape = NA, fill = c("brown1", "darkorange1")) +
  scale_x_discrete(labels = c("Relapse", "Remission")) +
  geom_jitter(width = 0.1, size = 0.7)+ 
  theme_classic() +
  labs(x = element_blank(), y = expression(paste("GFR (mL/min/1.73m" ^ "2", ")")))
dev.off()

#Stats between PR3 and Time
pr3$timepoint <- as.character(pr3$timepoint)
pr3_t <- pr3 %>% 
  t_test(pr3_titre ~ timepoint) %>% 
  add_significance()
pr3_t

pdf(file = "Blood - PR3 Titre over time.pdf", height = 3, width = 4)
ggplot(pr3, aes(x = timepoint, y = pr3_titre)) +
  geom_boxplot(outlier.shape = NA, fill = c("brown1", "darkorange1")) +
  scale_x_discrete(labels = c("Relapse", "Remission")) +
  geom_jitter(width = 0.1, size = 0.7)+ 
  theme_classic() +
  labs(x = element_blank(), y = "PR3 Titre (iu/mL)") +
  stat_pvalue_manual(pr3_t, label = "p.signif", tip.length = 0.01, y.position = 200)
dev.off()

#BVAS and CRP correlation
bvas_crp <- subset(crp_magenta, !is.na(bvas)) 

pdf(file = "Blood - CRP and BVAS correlation.pdf", height = 3.5, width = 5)
ggplot(bvas_crp, aes(x = bvas_score, y = crp)) +
  geom_point(aes(color = timepoint, shape = timepoint, fill = timepoint),
             color = "black", size = 2) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = "Timepoint (months)", 
                     labels = c("Zero", "Three")) +
  scale_fill_manual(values = c("brown1", "darkorange1"),
                    name = "Timepoint (months)", 
                    labels = c("Zero", "Three")) +
  geom_smooth(method = lm, se = FALSE, colour = "skyblue1") +
  labs(x = "BVAS Score", y = "CRP (mg/dL)") +
  theme_classic() +
  theme(legend.position = "top") +
  stat_cor(method = "pearson", label.x = 12, label.y = 60)
dev.off()
