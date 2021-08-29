##Differential expression between timepoints

library(limma)

load(file = "Vascmicroarray-final.RData")

groups = pData(vascFinalSet)$timepoint_short

f = factor(groups, levels = 
             c("untimed", "zero", "three", "twelve"))
design = model.matrix(~ 0 + f)
colnames(design) = c("untimed", "zero", "three", "twelve")
data.fit = lmFit(vascFinalSet, design)
contrast.matrix = makeContrasts(untimed-zero,
                                untimed-three,
                                untimed-twelve,
                                zero-three,
                                zero-twelve,
                                three-twelve,
                                levels = design)
data.fit.con = contrasts.fit(data.fit, contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

data.fit.eb$coefficients[1:10]

deg_timepoint <- topTable(data.fit.eb, number = Inf)
head(deg_timepoint)
hist(deg_timepoint$P.Value, col = "red",
     main = "Timepoints", xlab = "p-values")

##Differential expression - active_vs_hc

#Create data set with only untimed and zero samples
active_vs_hc_samples <- colnames(vascFinalSet)[vascFinalSet@phenoData@data$"timepoint_short"=="untimed" |
                                                  vascFinalSet@phenoData@data$"timepoint_short"=="zero"]
active_vs_hc_filtered <- vascFinalSet[ ,active_vs_hc_samples]

groups = pData(active_vs_hc_filtered)$timepoint_short

f = factor(groups, levels = 
             c("untimed", "zero"))
design = model.matrix(~ 0 + f)
colnames(design) = c("untimed", "zero")
data.fit = lmFit(active_vs_hc_filtered, design)
contrast.matrix = makeContrasts(zero-untimed,
                                levels = design)
data.fit.con = contrasts.fit(data.fit, contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

data.fit.eb$coefficients[1:10]

results <- decideTests(data.fit.eb)
summary(results)

active_vs_hc_deg <- topTable(data.fit.eb, adjust="BH", number = Inf)
head(active_vs_hc_deg)
hist(active_vs_hc_deg$P.Value, col = "red",
     main = "Active vs HC", xlab = "p-values")

write.csv(active_vs_hc_deg, file = "active_vs_hc trait/active_vs_hc_deg.csv")


##Differential expression - active_vs_remission

#Create data set with only zero and three samples
active_vs_remiss_samples <- colnames(vascFinalSet)[vascFinalSet@phenoData@data$"timepoint_short"=="zero" |
                                                 vascFinalSet@phenoData@data$"timepoint_short"=="three"]
active_vs_remiss_filtered <- vascFinalSet[ ,active_vs_remiss_samples]

groups = pData(active_vs_remiss_filtered)$timepoint_short

f = factor(groups, levels = 
             c("zero", "three"))
design = model.matrix(~ 0 + f)
colnames(design) = c("zero", "three")
data.fit = lmFit(active_vs_remiss_filtered, design)
contrast.matrix = makeContrasts(three-zero,
                                levels = design)
data.fit.con = contrasts.fit(data.fit, contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

data.fit.eb$coefficients[1:10]

active_vs_remiss_deg <- topTable(data.fit.eb, number = Inf)
head(active_vs_remiss_deg)
hist(active_vs_remiss_deg$P.Value, col = "red",
     main = "Active vs Remission", xlab = "p-values")

write.csv(active_vs_remiss_deg, file = "active_vs_remission trait/active_vs_remission_deg.csv")


##Differential expression - active_vs_twelve

#Create data set with only zero and twelve samples
active_vs_twelve_samples <- colnames(vascFinalSet)[vascFinalSet@phenoData@data$"timepoint_short"=="zero" |
                                                     vascFinalSet@phenoData@data$"timepoint_short"=="twelve"]
active_vs_twelve_filtered <- vascFinalSet[ ,active_vs_twelve_samples]

groups = pData(active_vs_twelve_filtered)$timepoint_short

f = factor(groups, levels = 
             c("zero", "twelve"))
design = model.matrix(~ 0 + f)
colnames(design) = c("zero", "twelve")
data.fit = lmFit(active_vs_twelve_filtered, design)
contrast.matrix = makeContrasts(twelve-zero,
                                levels = design)
data.fit.con = contrasts.fit(data.fit, contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

data.fit.eb$coefficients[1:10]

active_vs_twelve_deg <- topTable(data.fit.eb, number = Inf)
head(active_vs_twelve_deg)
hist(active_vs_twelve_deg$P.Value, col = "red",
     main = "Active vs Twelve", xlab = "p-values")

write.csv(active_vs_twelve_deg, file = "active_vs_twelve trait/active_vs_twelve_deg.csv")
