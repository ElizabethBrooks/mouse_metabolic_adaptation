#!/usr/bin/env Rscript

# R script to perform variance partitioning

# turn off scientific notation
options(scipen = 999)

# import libraries
library(edgeR)
library(tibble)
library(dplyr)
library(variancePartition)
library(BiocParallel)

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DA_12Nov2026/WAT_variance_partitioning_domesticus"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/MackLab/metabolic_adaptation/counted/counts_merged_samples.csv", row.names="gene")

# import grouping factor
factors <- read.csv(file="/Users/bamflappy/MackLab/input_data/experimental_design.csv", row.names="sample")

# set input FDR cutoff
cutFDR <- 0.05

# set LFC cut
#cutFC <- 1.2
#cutLFC <- log2(1.2)
cutLFC <- 0

# subset to tissue specific data
factors <- factors[factors$tissue == "WAT",]
keep_cols <- rownames(factors)
inputTable <- inputTable %>% select(all_of(keep_cols))

# trim the data table to remove lines with counting statistics (htseq)
removeList <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
countsTable <- inputTable[!row.names(inputTable) %in% removeList,]

# remove castaneus samples
remove <- rownames(factors[grepl("CAST", factors$mouseline),])
countsTable <- select(all_of(countsTable), -c(remove))
factors <- factors[!grepl("CAST", factors$mouseline),]

# convert the grouping data into factors 
targets <- as.data.frame(lapply(factors, as.factor))
rownames(targets) <- rownames(factors)

# specify the reference level
targets$treatment <- relevel(targets$treatment, ref = "Ad_lib")
targets$mouseline <- relevel(targets$mouseline, ref = "SARA")
#targets$tissue <- relevel(targets$tissue, ref = "BAT")
targets$sex <- relevel(targets$sex, ref = "Male")

# sort samples alphabetically
countsTable <- countsTable[, order(colnames(countsTable))]
targets <- targets[order(rownames(targets)), ]

# check sample naming
setdiff(colnames(countsTable), rownames(targets))
setdiff(rownames(targets), colnames(countsTable))

# there are at least 1 samples per individual
smallestGroupSize <- 1

# keep only rows that have a count of at least 10 for a minimal number of samples
keep <- rowSums(cpm(countsTable) >= 10) >= smallestGroupSize
list <- DGEList(countsTable[keep, ])

# percentage of genes lost from filtering
(nrow(countsTable) - nrow(list)) / nrow(countsTable) * 100

# use TMM normalization to eliminate composition biases between libraries
#dge <- calcNormFactors(list)
dge <- normLibSizes(list)

# specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ treatment + sex + mouseline + treatment:mouseline + mouseline:sex

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ treatment + sex + mouseline + treatment:mouseline + mouseline:sex
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart_treatment_sex_mouseline_treatment.mouseline_mouseline.sex.jpg")
plotVarPart(sortCols(vp))
dev.off()

# The variable to be tested must be a fixed effect
form <- ~ treatment + sex + mouseline + treatment:mouseline

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ treatment + sex + mouseline + treatment:mouseline
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart_treatment_sex_mouseline_treatment.mouseline.jpg")
plotVarPart(sortCols(vp))
dev.off()
