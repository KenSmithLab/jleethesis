##Module Preservation

library(WGCNA)
library(ggplot2)
setwd("~/Jacinta/Blood/Module Preservation/")

options(stringsAsFactors = FALSE);

load("~/Jacinta/microMatrix.RData")
coloursMicro <- moduleColors

#Module Preservation for all Myeloid Cells
load("~/Jacinta/Blood/Pre_processing/myeloid_final - Blood.RData")
myeloid_matrix <- 
  t(as.matrix(GetAssayData(myeloid.final)))

gsg = goodSamplesGenes(myeloid_matrix, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(myeloid_matrix)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(myeloid_matrix)[!gsg$goodSamples], collapse = ", ")))
  myeloid_clean = myeloid_matrix[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = myeloid_clean));
multiColor = list(Micro = coloursMicro);

system.time( {
  total_mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 150,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(total_mp, file = "TotalCells_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(total_mp$quality$observed[[ref]][[test]][, -1], total_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(total_mp$quality$Z[[ref]][[test]][, -1], total_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(total_mp$preservation$observed[[ref]][[test]])
moduleSizes = total_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(total_mp$preservation$observed[[ref]][[test]][, 2], total_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Total_myeloid_cells_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();

##Module Preservation for cells of myeloid subsets
load("Myeloid_subsets_cells_matrix.RData")

#Classical
gsg = goodSamplesGenes(classical_matrix, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_matrix)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_matrix)[!gsg$goodSamples], collapse = ", ")))
  classical_cells_clean = classical_matrix[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_cells_clean));
multiColor = list(Micro = coloursMicro);

system.time( {
  classical_cells_mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 100,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(classical_cells_mp, file = "Classical_cells_modulePreservation.RData")


#Intermediate
gsg = goodSamplesGenes(intermediate_matrix, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(intermediate_matrix)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(intermediate_matrix)[!gsg$goodSamples], collapse = ", ")))
  intermediate_matrix_clean = intermediate_matrix[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = intermediate_matrix_clean));
multiColor = list(Micro = coloursMicro);

system.time( {
  intermediate_cells_mp = modulePreservation(multiExpr, multiColor,
                                          referenceNetworks = 1,
                                          nPermutations = 100,
                                          randomSeed = 1,
                                          quickCor = 0,
                                          verbose = 3)
})

save(intermediate_cells_mp, file = "Intermediate_cells_modulePreservation.RData")

#Non-Classical
gsg = goodSamplesGenes(nonclassical_matrix, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(nonclassical_matrix)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(nonclassical_matrix)[!gsg$goodSamples], collapse = ", ")))
  nonclassical_cells_clean = nonclassical_matrix[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = nonclassical_cells_clean));
multiColor = list(Micro = coloursMicro);

system.time( {
  nonclassical_cells_mp = modulePreservation(multiExpr, multiColor,
                                             referenceNetworks = 1,
                                             nPermutations = 100,
                                             randomSeed = 1,
                                             quickCor = 0,
                                             verbose = 3)
})

save(nonclassical_cells_mp, file = "Nonclassical_cells_modulePreservation.RData")

#Myeloid DC
gsg = goodSamplesGenes(myeloid_dc_matrix, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(myeloid_dc_matrix)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(myeloid_dc_matrix)[!gsg$goodSamples], collapse = ", ")))
  myeloid_dc_cells_clean = myeloid_dc_matrix[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = myeloid_dc_cells_clean));
multiColor = list(Micro = coloursMicro);

system.time( {
  myeloid_dc_cells_mp = modulePreservation(multiExpr, multiColor,
                                             referenceNetworks = 1,
                                             nPermutations = 100,
                                             randomSeed = 1,
                                             quickCor = 0,
                                             verbose = 3)
})

save(myeloid_dc_cells_mp, file = "MyeloidDC_cells_modulePreservation.RData")












#Module Preservation for all average of Myeloid cells per patient
load("myeloid_all_matrix.RData")
myeloid_all_matrix <- t(myeloid_all_matrix)

#Check and isolate genes with expression
gsg = goodSamplesGenes(myeloid_all_matrix, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(myeloid_all_matrix)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(myeloid_all_matrix)[!gsg$goodSamples], collapse = ", ")))
  myeloid_all_clean = myeloid_all_matrix[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = myeloid_all_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(mp, file = "TotalMyeloid_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Total_myeloid_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();

pdf(file = "Total_myeloid_mp_plot.pdf", width = 5, height = 6)
plot(moduleSizes[plotMods], plotData[plotMods, 1], col = 1, bg = modColors[plotMods], pch = 21,
     main = mains[1],
     cex = 2.4,
     ylab = mains[1], xlab = "Module size", log = "x",
     ylim = ylim,
     xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.5);




#Module Preservation of individual myeloid subsets that have been averaged per sample
load("Myeloid_subsets_matrix.RData")

#Classical monocyte subset
classical <- t(classical)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical)[!gsg$goodSamples], collapse = ", ")))
  classical_clean = classical[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(classical_mp, file = "Classical_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_mp$quality$observed[[ref]][[test]][, -1], classical_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_mp$quality$Z[[ref]][[test]][, -1], classical_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_mp$preservation$observed[[ref]][[test]][, 2], classical_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Classical_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


##Intermediate monocyte subset
intermediate <- t(intermediate)

#Check and isolate genes with expression
gsg = goodSamplesGenes(intermediate, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(intermediate)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(intermediate)[!gsg$goodSamples], collapse = ", ")))
  intermediate_clean = intermediate[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = intermediate_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  intermediate_mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(intermediate_mp, file = "Intermediate_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(intermediate_mp$quality$observed[[ref]][[test]][, -1], intermediate_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(intermediate_mp$quality$Z[[ref]][[test]][, -1], intermediate_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(intermediate_mp$preservation$observed[[ref]][[test]])
moduleSizes = intermediate_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(intermediate_mp$preservation$observed[[ref]][[test]][, 2], intermediate_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Intermediate_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


#Non-Classical monocyte subset
nonclassical <- t(nonclassical)

#Check and isolate genes with expression
gsg = goodSamplesGenes(nonclassical, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(nonclassical)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(nonclassical)[!gsg$goodSamples], collapse = ", ")))
  nonclassical_clean = nonclassical[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = nonclassical_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  nonclassical_mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(nonclassical_mp, file = "NonClassical_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(nonclassical_mp$quality$observed[[ref]][[test]][, -1], nonclassical_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(nonclassical_mp$quality$Z[[ref]][[test]][, -1], nonclassical_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(nonclassical_mp$preservation$observed[[ref]][[test]])
moduleSizes = nonclassical_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(nonclassical_mp$preservation$observed[[ref]][[test]][, 2], nonclassical_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Nonclassical_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


#MyeloidDC monocyte subset
myeloid_dc <- t(myeloid_dc)

#Check and isolate genes with expression
gsg = goodSamplesGenes(myeloid_dc, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(myeloid_dc)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(myeloid_dc)[!gsg$goodSamples], collapse = ", ")))
  myeloid_dc_clean = myeloid_dc[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = myeloid_dc_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  myeloid_dc_mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(myeloid_dc_mp, file = "MyeloidDC_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(myeloid_dc_mp$quality$observed[[ref]][[test]][, -1], myeloid_dc_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(myeloid_dc_mp$quality$Z[[ref]][[test]][, -1], myeloid_dc_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(myeloid_dc_mp$preservation$observed[[ref]][[test]])
moduleSizes = myeloid_dc_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(myeloid_dc_mp$preservation$observed[[ref]][[test]][, 2], myeloid_dc_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Myeloid_dc_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();

#Remove modules that are not in all subsets, dark turquoise and lightyellow
load("Classical_modulePreservation.RData")
load("Intermediate_modulePreservation.RData")
load("NonClassical_modulePreservation.RData")
load("MyeloidDC_modulePreservation.RData")
ref = 1
test = 2
modColors = rownames(classical_mp$preservation$observed[[ref]][[test]])

Z_classical <- classical_mp$preservation$Z[[ref]][[test]][, 2]
Z_intermediate <- intermediate_mp$preservation$Z[[ref]][[test]][, 2]
Z_nonclassical <- nonclassical_mp$preservation$Z[[ref]][[test]][, 2]
Z_myeloid_dc <- myeloid_dc_mp$preservation$Z[[ref]][[test]][, 2]
Z_sum <- as.data.frame(cbind(modColors, Z_classical, Z_intermediate, Z_nonclassical, Z_myeloid_dc))

Z_long <- Z_sum %>% 
  pivot_longer(cols = Z_classical:Z_myeloid_dc,
               names_to = "cluster",
               values_to = "Zsum")

#Get magenta Z summary values for each cluster
Z_magenta <- Z_long[grep("magenta", Z_long$modColors), ]
Z_magenta$cluster <- factor(Z_magenta$cluster, 
                            levels = c("Z_classical", "Z_nonclassical", "Z_intermediate", "Z_myeloid_dc"))
Z_magenta$Zsum <- as.numeric(Z_magenta$Zsum)

pdf(file = "Magenta Zsummary for myeloid subsets.pdf", width = 5, height = 3)
ggplot(Z_magenta, aes(x = cluster, y = Zsum, fill = cluster)) +
  geom_col() +
  scale_x_discrete(labels = c("Classical","Non-Classical", "Intermediate", "Myeloid DC")) +
  theme_bw() +
  labs(x = "Subset", y = "Magenta Preservation Z Summary") +
  theme(legend.position = "none") +
  geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

#Get pink Z summary values for each cluster
Z_pink <- Z_long[grep("pink", Z_long$modColors), ]
Z_pink$cluster <- factor(Z_pink$cluster, 
                         levels = c("Z_classical", "Z_nonclassical", "Z_intermediate", "Z_myeloid_dc"))
Z_pink$Zsum <- as.numeric(Z_pink$Zsum)

pdf(file = "Pink Zsummary for myeloid subsets.pdf", width = 5, height = 3)
ggplot(Z_pink, aes(x = cluster, y = Zsum, fill = cluster)) +
  geom_col() +
  scale_x_discrete(labels = c("Classical","Non-Classical", "Intermediate", "Myeloid DC")) +
  theme_bw() +
  labs(x = "Subset", y = "Pink Preservation Z Summary") +
  theme(legend.position = "none") +
  geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()




#Get a module score for each cell type for the Magenta Module
magentaGenes <- as.vector(read.csv("~/Jacinta/magentaGenes.csv"))
magentaGenes <- magentaGenes$x
load("Myeloid_subsets_seurat.RData")

classical_cells <- AddModuleScore(classical_cells, list(magentaGenes), name = "magentaGenes")
intermediate_cells <- AddModuleScore(intermediate_cells, list(magentaGenes), name = "magentaGenes")
nonclassical_cells <- AddModuleScore(nonclassical_cells, list(magentaGenes), name = "magentaGenes")
myeloid_dc_cells <- AddModuleScore(myeloid_dc_cells, list(magentaGenes), name = "magentaGenes")
myeloid.final <- AddModuleScore(myeloid.final, list(magentaGenes), name = "magentaGenes")

FeaturePlot(myeloid.final, features = "magentaGenes1", order = TRUE)

p1 <- FeaturePlot(classical_cells, features = "magentaGenes1", order = TRUE)
p2 <- FeaturePlot(intermediate_cells, features = "magentaGenes1", order = TRUE)
p3 <- FeaturePlot(nonclassical_cells, features = "magentaGenes1", order = TRUE)
p4 <- FeaturePlot(myeloid_dc_cells, features = "magentaGenes1", order = TRUE)
p5 <- p1|p2|p3|p4
p5





#Module Preservation for average of Classical cell subsets per patient
load("~/Jacinta/Blood/Pre_processing/Classical_subsets_avg_matrix.RData")

#Classical 0
classical_0 <- t(classical_0)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_0)[!gsg$goodSamples], collapse = ", ")))
  classical_0_clean = classical_0[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_0_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_0_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_0_mp, file = "classical_0_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_0_mp$quality$observed[[ref]][[test]][, -1], classical_0_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_0_mp$quality$Z[[ref]][[test]][, -1], classical_0_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_0_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_0_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_0_mp$preservation$observed[[ref]][[test]][, 2], classical_0_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_0_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();

#Classical 1
classical_1 <- t(classical_1)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_1, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_1)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_1)[!gsg$goodSamples], collapse = ", ")))
  classical_1_clean = classical_1[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_1_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_1_mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
})

save(classical_1_mp, file = "classical_1_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_1_mp$quality$observed[[ref]][[test]][, -1], classical_1_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_1_mp$quality$Z[[ref]][[test]][, -1], classical_1_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_1_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_1_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_1_mp$preservation$observed[[ref]][[test]][, 2], classical_1_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_1_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


#Classical 2
classical_2 <- t(classical_2)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_2, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_2)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_2)[!gsg$goodSamples], collapse = ", ")))
  classical_2_clean = classical_2[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_2_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_2_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_2_mp, file = "classical_2_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_2_mp$quality$observed[[ref]][[test]][, -1], classical_2_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_2_mp$quality$Z[[ref]][[test]][, -1], classical_2_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_2_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_2_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_2_mp$preservation$observed[[ref]][[test]][, 2], classical_2_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_2_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();



#Classical 3
classical_3 <- t(classical_3)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_3, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_3)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_3)[!gsg$goodSamples], collapse = ", ")))
  classical_3_clean = classical_3[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_3_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_3_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_3_mp, file = "classical_3_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_3_mp$quality$observed[[ref]][[test]][, -1], classical_3_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_3_mp$quality$Z[[ref]][[test]][, -1], classical_3_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_3_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_3_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_3_mp$preservation$observed[[ref]][[test]][, 2], classical_3_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_3_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();



#Classical 4
classical_4 <- t(classical_4)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_4, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_4)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_4)[!gsg$goodSamples], collapse = ", ")))
  classical_4_clean = classical_4[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_4_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_4_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_4_mp, file = "classical_4_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_4_mp$quality$observed[[ref]][[test]][, -1], classical_4_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_4_mp$quality$Z[[ref]][[test]][, -1], classical_4_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_4_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_4_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_4_mp$preservation$observed[[ref]][[test]][, 2], classical_4_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_4_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();



#Classical 7
classical_7 <- t(classical_7)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_7, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_7)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_7)[!gsg$goodSamples], collapse = ", ")))
  classical_7_clean = classical_7[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_7_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_7_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_7_mp, file = "classical_7_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_7_mp$quality$observed[[ref]][[test]][, -1], classical_7_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_7_mp$quality$Z[[ref]][[test]][, -1], classical_7_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_7_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_7_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_7_mp$preservation$observed[[ref]][[test]][, 2], classical_7_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_7_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


#Classical 9
classical_8 <- t(classical_8)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_8, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_8)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_8)[!gsg$goodSamples], collapse = ", ")))
  classical_8_clean = classical_8[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_8_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_8_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_8_mp, file = "classical_8_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_8_mp$quality$observed[[ref]][[test]][, -1], classical_8_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_8_mp$quality$Z[[ref]][[test]][, -1], classical_8_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_8_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_8_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_8_mp$preservation$observed[[ref]][[test]][, 2], classical_8_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_8_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();



#Classical 10
classical_10 <- t(classical_10)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_10, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_10)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_10)[!gsg$goodSamples], collapse = ", ")))
  classical_10_clean = classical_10[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_10_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_10_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_10_mp, file = "classical_10_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_10_mp$quality$observed[[ref]][[test]][, -1], classical_10_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_10_mp$quality$Z[[ref]][[test]][, -1], classical_10_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_10_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_10_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_10_mp$preservation$observed[[ref]][[test]][, 2], classical_10_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_10_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


#Classical 11
classical_11 <- t(classical_11)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_11, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_11)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_11)[!gsg$goodSamples], collapse = ", ")))
  classical_11_clean = classical_11[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_11_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_11_mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
})

save(classical_11_mp, file = "classical_11_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_11_mp$quality$observed[[ref]][[test]][, -1], classical_11_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_11_mp$quality$Z[[ref]][[test]][, -1], classical_11_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_11_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_11_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_11_mp$preservation$observed[[ref]][[test]][, 2], classical_11_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_11_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


#Classical 13
classical_13 <- t(classical_13)

#Check and isolate genes with expression
gsg = goodSamplesGenes(classical_13, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(classical_13)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(classical_13)[!gsg$goodSamples], collapse = ", ")))
  classical_13_clean = classical_13[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = classical_13_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  classical_13_mp = modulePreservation(multiExpr, multiColor,
                                       referenceNetworks = 1,
                                       nPermutations = 200,
                                       randomSeed = 1,
                                       quickCor = 0,
                                       verbose = 3)
})

save(classical_13_mp, file = "classical_13_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(classical_13_mp$quality$observed[[ref]][[test]][, -1], classical_13_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(classical_13_mp$quality$Z[[ref]][[test]][, -1], classical_13_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(classical_13_mp$preservation$observed[[ref]][[test]])
moduleSizes = classical_13_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(classical_13_mp$preservation$observed[[ref]][[test]][, 2], classical_13_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="classical_13_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();

#Remove modules that are not in all subsets, dark turquoise and lightyellow
load("classical_0_modulePreservation.RData")
load("classical_1_modulePreservation.RData")
load("classical_2_modulePreservation.RData")
load("classical_3_modulePreservation.RData")
load("classical_4_modulePreservation.RData")
load("classical_7_modulePreservation.RData")
load("classical_8_modulePreservation.RData")
load("classical_10_modulePreservation.RData")
load("classical_11_modulePreservation.RData")
load("classical_13_modulePreservation.RData")
ref = 1
test = 2
modColors13 = rownames(classical_13_mp$preservation$observed[[ref]][[test]])


Z_0 <- classical_0_mp$preservation$Z[[ref]][[test]][, 2]
Z_0 <- Z_0[-c(7,15)]
Z_1 <- classical_1_mp$preservation$Z[[ref]][[test]][, 2]
Z_1 <- Z_1[-c(7,15)]
Z_2 <- classical_2_mp$preservation$Z[[ref]][[test]][, 2]
Z_2 <- Z_2[-c(7,15)]
Z_3 <- classical_3_mp$preservation$Z[[ref]][[test]][, 2]
Z_3 <- Z_3[-c(7,15)]
Z_4 <- classical_4_mp$preservation$Z[[ref]][[test]][, 2]
Z_4 <- Z_4[-c(7,15)]
Z_7 <- classical_7_mp$preservation$Z[[ref]][[test]][, 2]
Z_7 <- Z_7[-c(7,15)]
Z_8 <- classical_8_mp$preservation$Z[[ref]][[test]][, 2]
Z_8 <- Z_8[-c(7,15)]
Z_10 <- classical_10_mp$preservation$Z[[ref]][[test]][, 2]
Z_10 <- Z_10[-c(7,15)]
Z_11 <- classical_11_mp$preservation$Z[[ref]][[test]][, 2]
Z_11 <- Z_11[-c(7,15)]
Z_13 <- classical_13_mp$preservation$Z[[ref]][[test]][, 2]
Z_sum <- as.data.frame(cbind(modColors13, Z_0, Z_1, Z_2, Z_3, Z_4, 
                             Z_7, Z_8, Z_10, Z_11, Z_13))

library(tidyr)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(13)
subset_cols = cols[c(1, 2, 3, 4, 5, 6, 9, 10, 12, 13, 15)]

Z_long <- Z_sum %>% 
  pivot_longer(cols = Z_0:Z_13,
               names_to = "cluster",
               values_to = "Zsum")

#Get magenta Z summary values for each cluster
Z_magenta <- Z_long[grep("magenta", Z_long$modColors13), ]
Z_magenta$cluster <- factor(Z_magenta$cluster, 
                               levels = c("Z_0", "Z_1", "Z_2", "Z_3", "Z_4", "Z_7", 
                                          "Z_8", "Z_10", "Z_11", "Z_13"))
Z_magenta$Zsum <- as.numeric(Z_magenta$Zsum)

pdf(file = "Magenta Zsummary for classical subset.pdf", width = 5, height = 3)
ggplot(Z_magenta, aes(x = cluster, y = Zsum, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = subset_cols) +
  theme_bw() +
  labs(x = "Cluster", y = "Magenta Preservation Z Summary") +
  theme(legend.position = "none") +
  geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

#Get pink Z summary values for each cluster
Z_pink <- Z_long[grep("pink", Z_long$modColors13), ]
Z_pink$cluster <- factor(Z_pink$cluster, 
                            levels = c("Z_0", "Z_1", "Z_2", "Z_3", "Z_4", "Z_7", 
                                       "Z_8", "Z_10", "Z_11", "Z_13"))
Z_pink$Zsum <- as.numeric(Z_pink$Zsum)

pdf(file = "Pink Zsummary for classical subset.pdf", width = 5, height = 3)
ggplot(Z_pink, aes(x = cluster, y = Zsum, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = subset_cols) +
  theme_bw() +
  labs(x = "Cluster", y = "Pink Preservation Z Summary") +
  theme(legend.position = "none") +
  geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()





#Module Preservation of timepoints that have been averaged per sample
load("Myeloid_timepoints_matrix.RData")

#Control
control <- t(control)

#Check and isolate genes with expression
gsg = goodSamplesGenes(control, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(control)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(control)[!gsg$goodSamples], collapse = ", ")))
  control_clean = control[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = control_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  control_mp = modulePreservation(multiExpr, multiColor,
                                  referenceNetworks = 1,
                                  nPermutations = 200,
                                  randomSeed = 1,
                                  quickCor = 0,
                                  verbose = 3)
})

save(control_mp, file = "control_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(control_mp$quality$observed[[ref]][[test]][, -1], control_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(control_mp$quality$Z[[ref]][[test]][, -1], control_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(control_mp$preservation$observed[[ref]][[test]])
moduleSizes = control_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(control_mp$preservation$observed[[ref]][[test]][, 2], control_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="control_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


##Zero timepoint
zero <- t(zero)

#Check and isolate genes with expression
gsg = goodSamplesGenes(zero, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(zero)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(zero)[!gsg$goodSamples], collapse = ", ")))
  zero_clean = zero[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = zero_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  zero_mp = modulePreservation(multiExpr, multiColor,
                               referenceNetworks = 1,
                               nPermutations = 200,
                               randomSeed = 1,
                               quickCor = 0,
                               verbose = 3)
})

save(zero_mp, file = "zero_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(zero_mp$quality$observed[[ref]][[test]][, -1], zero_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(zero_mp$quality$Z[[ref]][[test]][, -1], zero_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(zero_mp$preservation$observed[[ref]][[test]])
moduleSizes = zero_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(zero_mp$preservation$observed[[ref]][[test]][, 2], zero_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="zero_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();


#Three month subset
three <- t(three)

#Check and isolate genes with expression
gsg = goodSamplesGenes(three, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(three)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(three)[!gsg$goodSamples], collapse = ", ")))
  three_clean = three[gsg$goodSamples, gsg$goodGenes]
}

setLabels = c("Micro", "Single")
multiExpr = list(Micro = list(data = microMatrix), Single = list(data = three_clean));
multiColor = list(Micro = coloursMicro);

#Run module preservation function
system.time( {
  three_mp = modulePreservation(multiExpr, multiColor,
                                referenceNetworks = 1,
                                nPermutations = 200,
                                randomSeed = 1,
                                quickCor = 0,
                                verbose = 3)
})

save(three_mp, file = "three_modulePreservation.RData")

#Analysis of module preservation results
#Isolate statistics and Z scores
ref = 1
test = 2
statsObs = cbind(three_mp$quality$observed[[ref]][[test]][, -1], three_mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(three_mp$quality$Z[[ref]][[test]][, -1], three_mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#Plot the preservation medianRank and Zsummary
# Module labels and module sizes are also contained in the results
modColors = rownames(three_mp$preservation$observed[[ref]][[test]])
moduleSizes = three_mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(three_mp$preservation$observed[[ref]][[test]][, 2], three_mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="three_avg_mp_plot.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();



#Visualise Z summary for as a barplot
load("control_modulePreservation.RData")
load("zero_modulePreservation.RData")
load("three_modulePreservation.RData")
ref = 1
test = 2
modColors = rownames(control_mp$preservation$observed[[ref]][[test]])

Z_control <- control_mp$preservation$Z[[ref]][[test]][, 2]
Z_zero <- zero_mp$preservation$Z[[ref]][[test]][, 2]
Z_three <- three_mp$preservation$Z[[ref]][[test]][, 2]
Z_sum <- as.data.frame(cbind(modColors, Z_control, Z_zero, Z_three))

Z_long <- Z_sum %>% 
  pivot_longer(cols = Z_control:Z_three,
               names_to = "cluster",
               values_to = "Zsum")

#Get magenta Z summary values for each cluster
Z_magenta <- Z_long[grep("magenta", Z_long$modColors), ]
Z_magenta$cluster <- factor(Z_magenta$cluster, 
                            levels = c("Z_control", "Z_zero", "Z_three"))
Z_magenta$Zsum <- as.numeric(Z_magenta$Zsum)

pdf(file = "Blood - Magenta Zsummary over time.pdf", width = 5, height = 3)
ggplot(Z_magenta, aes(x = cluster, y = Zsum, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = c("skyblue1","brown1", "darkorange1")) +
  scale_x_discrete(labels = c("HC", "Relapse", "Remission")) +
  theme_bw() +
  labs(x = element_blank(), y = "Magenta Preservation Z Summary") +
  theme(legend.position = "none") +
  geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()

#Get pink Z summary values for each cluster
Z_pink <- Z_long[grep("pink", Z_long$modColors), ]
Z_pink$cluster <- factor(Z_pink$cluster, 
                         levels = c("Z_control", "Z_zero", "Z_three"))
Z_pink$Zsum <- as.numeric(Z_pink$Zsum)

pdf(file = "Blood - Pink Zsummary over time.pdf", width = 5, height = 3)
ggplot(Z_pink, aes(x = cluster, y = Zsum, fill = cluster)) +
  geom_col() +
  scale_fill_manual(values = c("skyblue1","brown1", "darkorange1")) +
  scale_x_discrete(labels = c("HC", "Relapse", "Remission")) +
  theme_bw() +
  labs(x = element_blank(), y = "Pink Preservation Z Summary") +
  theme(legend.position = "none") +
  geom_hline(yintercept=2, linetype="dashed", color = "red")
dev.off()







