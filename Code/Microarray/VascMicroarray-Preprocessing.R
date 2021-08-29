##Vasculitis microarray data pre processing

library(arrayQualityMetrics)
library(oligo)
library(limma)
library(gcrma)
library(sva)
library(WGCNA)
library(hugene11stprobeset.db)
library(hugene11sttranscriptcluster.db)
library(devtools)
library(dplyr)
library(Biobase)
library(ggfortify)

setwd("~/VascMicroarray/Usable Data")

list = list.files(full.names = TRUE, recursive = T, patt ="*.CEL")

#Checking samples for adequacy
ctrllist <- list.files(full.names = TRUE, recursive = TRUE, patt ="*.CEL")
ctrldata <- read.celfiles(ctrllist)
summary(exprs(ctrldata))
boxplot(ctrldata, target = "probeset")
ctrlrma <- rma(ctrldata)

vasclist <- list.files(full.names = TRUE, recursive = TRUE, patt ="*.CEL")
vascdata <- read.celfiles(vasclist)
boxplot(vascdata, target = "probeset")
summary(exprs(vascdata))
vascrma <- rma(vascdata)


#Read in phenotype data
pheno <- read.csv("Phenodata.csv")
row.names(pheno) = pheno$filename
pheno <- AnnotatedDataFrame(pheno)

rawdata <- read.celfiles(list, phenoData = pheno)


#create expression set
expset <- oligo::rma(rawdata)

#create data frame with just expression values
expmatrix <- exprs(expset) 

#Microarray Pictures
ph <- rawdata@phenoData
for (i in 1:6)
{
  name = paste("image",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=ph@data$sample[i])
  dev.off()
}

#Histogram
for (i in 1:6)
{
  name = paste("histogram",i,".jpg",sep="")
  jpeg(name)
  oligo::hist(rawdata[,i], target = "probeset", 
              lwd=2,which='pm',
              ylab='Density',
              xlab='Log2 intensities',main=ph@data$sample[i])
  dev.off()
}

#Boxplot of rawdata
boxplot(rawdata, target = "probeset")

#Boxplot of normalized data
boxplot(expset, target = "probeset")

#PCA of normalized data

PC_batch = prcomp(t(expmatrix),scale.=TRUE)
batch_name <- as.factor(pheno$batch_name)

PC_batch_plot <- as.data.frame(PC_batch$x)
PC_batch_plot <- cbind(PC_batch_plot, batch_name)
PC_batch_plot <- subset(PC_batch_plot, PC1 < 200)
ggplot(PC_batch_plot ,aes(x=PC1, y=PC2, color = batch_name )) +
  geom_point() +
  scale_color_discrete(name = "Batch") +
  theme_classic()
ggsave("PCplot.pdf", height = 4, width = 6)


#Array quality metrics
arrayQualityMetrics(expressionset = expset,
                    outdir = "Report_for_Vascmicroarray",
                    force = TRUE)
#Outliers that were high for 2 or more qualities was 
#deleted from original sample

#Batch correction is not needed
#Correcting for batch effects
#ph <- rawdata@phenoData
#batchcorrdata = ComBat(dat = expmatrix,
                       #batch = ph@data$batch_name,
                      # mod = NULL,
                      # par.prior = TRUE,
                       #prior.plots = FALSE)

#save(batchcorrdata, file = "batchcorrdata.Rdata")

#Samples not batch corrected as PC plot didn't show batch differentiation

#Filtering lowly expressed genes
vasc_medians <- rowMedians(expmatrix)

man_threshold <- 4.25

pdf(file = "medianHist.pdf", height = 8, width = 15)
hist_medians <- hist(vasc_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)
dev.off()

#Transcripts that do not have intensities larger than the threshold in at 
#least as many arrays as the smallest experimental group are excluded
no_of_samples <- table(paste0(pData(expset)$cohort))

samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(expmatrix, 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
                              table(idx_man_threshold)

vasc_manfiltered <- subset(expset, idx_man_threshold)
vasc_filtermatrix <- exprs(vasc_manfiltered)


#Transpose
trans_expmatrix <- t(vasc_filtermatrix)

#Find probes with matching gene annotation
library(hugene11sttranscriptcluster.db)

columns(hugene11sttranscriptcluster.db)

keys = keys(hugene11sttranscriptcluster.db)
annots <- AnnotationDbi::select(
  hugene11sttranscriptcluster.db,
  keys    = keys,
  columns = c("GENENAME", "SYMBOL", "ENTREZID"),
  keytype = "PROBEID"
)

head(unique(annots$GENENAME))
#Remove NA in annots
annots = subset(annots, !is.na(SYMBOL))
head(annots)

#Keep only 1 row if there are duplicates
annots <- as.data.frame(annots)
annots <- unique(annots)

#Delete rows that have multiple probe ID - 
#rows that have multiple genes for 1 probe
annots <- annots[!(duplicated(annots$PROBEID) | duplicated(annots$PROBEID, fromLast = TRUE)), ]

#Change row names of annots to probeID
rownames(annots) <- annots$PROBEID

#Remove genes with multiple probes using collapseRows and 
#selecting probes with the highest variance
cleanAnnots <- collapseRows(datET = vasc_filtermatrix,
                            rowGroup = annots$SYMBOL,
                            rowID = unique(annots$PROBEID),
                            method = "maxRowVariance")

#Extract a final matrix containing Gene and selected probeID
usefulProbes <- as.data.frame(cleanAnnots$group2row)
annotsMatches <- Reduce(intersect, list(usefulProbes$selectedRowID, rownames(annots)))
finalAnnots <- annots[annotsMatches, ]

write.csv(finalAnnots, "annots.csv")

#Match useful probes with expression set assay data
ids_to_include <- (featureNames(vasc_manfiltered) %in% usefulProbes$selectedRowID)
vascFinalSet <- subset(vasc_manfiltered, ids_to_include)
validObject(vascFinalSet)

#Match useful probes with expression set feature data
fData(vascFinalSet)$PROBEID <- rownames(fData(vascFinalSet))
fData(vascFinalSet) <- left_join(fData(vascFinalSet), finalAnnots)
rownames(fData(vascFinalSet)) <- fData(vascFinalSet)$PROBEID 
validObject(vascFinalSet)

vascFinalMatrix <- exprs(vascFinalSet)


#Match useful probes  with annotations with matrix
#trans_expmatrix = as.data.frame(trans_expmatrix)
#matches = Reduce(intersect, list(usefulProbes$selectedRowID, names(trans_expmatrix)))
#wgcnaCleanMatrix = trans_expmatrix[, matches]

#wgcnaCleanMatrix = wgcnaCleanMatrix[order(row.names(wgcnaCleanMatrix)), ]

#cleanMatrix = t(wgcnaCleanMatrix)

#Boxplot of normalized data
boxplot(vascFinalSet, target = "probeset")



save(vascFinalSet, vascFinalMatrix, finalAnnots, file = "Vascmicroarray-final.RData")

PC_batch = prcomp(t(vascFinalMatrix),scale.=TRUE)
batch_name <- as.factor(pheno$batch_name)
plot(PC_batch$x[ ,1:2], col = batch_name)
