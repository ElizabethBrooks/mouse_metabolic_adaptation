#!/usr/bin/env Rscript

# R script to perform DE analysis with DREAM

# turn off scientific notation
options(scipen = 999)

# import libraries
library(DESeq2)
library(tibble)
library(rcartocolor)
library(ggVennDiagram)
library(dplyr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(ggplotify)
library(ggpubr)
library(RColorBrewer)

##
# Inputs & Data Prep
##

# set working directory for outputs
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DE_19Nov2026/WAT_domesticus_LFC0"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/MackLab/metabolic_adaptation/counted/counts_merged_samples.csv", row.names="gene")

# import grouping factor
factors <- read.csv(file="/Users/bamflappy/MackLab/input_data/experimental_design.csv", row.names="sample")

# set input FDR cutoff
cutFDR <- 0.05

# set LFC cut
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

# verify that the order of the samples in the counts and groupings files match
#colnames(countsTable)
#rownames(factors)

# convert the grouping data into factors 
factors$individual <- factor(factors$individual)
factors$mouseline <- factor(factors$mouseline)
factors$sex <- factor(factors$sex)
factors$treatment <- factor(factors$treatment)
factors$tissue <- factor(factors$tissue)

# create DESeqDataSet list object
dds <- DESeqDataSetFromMatrix(countData = countsTable,
                              colData = factors,
                              design = ~ treatment + mouseline + sex + treatment:mouseline + sex:mouseline)

# inspect the list object
#dds

# specify the reference levels
dds$treatment <- relevel(dds$treatment, ref = "Ad_lib")
dds$mouseline <- relevel(dds$mouseline, ref = "MANF")
dds$sex <- relevel(dds$sex, ref = "Male")

##
# Pre-Filtering
##

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])
scales::show_col(plotColors)

# there are at least 1 samples per individual
smallestGroupSize <- 1

# keep only rows that have a count of at least 10 for a minimal number of samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

# filter the list object
dds <- dds[keep,]

# check out the number of kept genes
nrow(dds)

##
# Data Exploration
##

# vst the data
vsd <- vst(dds, blind=FALSE)

# visualize the overall effect of the experimental conditions or any batch effects
plotPCA(vsd, intgroup=c("treatment", "mouseline"))

# save the PCA
pcaData <- plotPCA(vsd, intgroup=c("treatment", "mouseline"), returnData=FALSE)

# customize the PCA
ggplot(pcaData@data, aes(PC1, PC2, color=sex, shape=mouseline, fill=treatment)) +
  geom_point(size=3, stroke=1.25) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() + 
  scale_shape_manual(values = c(21, 24)) +
  scale_colour_manual(values = c(plotColors[4], plotColors[3])) +
  scale_fill_manual(values = c(plotColors[1], plotColors[2])) +
  guides(color=guide_legend(override.aes=list(shape=21))) +
  guides(fill=guide_legend(override.aes=list(shape=21,fill=c(plotColors[1], plotColors[2]),color=c("white","white")))) +
  theme_minimal()

# store the PCA plot
sample_pca <- ggplot(pcaData@data, aes(PC1, PC2, color=sex, shape=mouseline, fill=treatment)) +
  geom_point(size=3, stroke=1.25) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() + 
  scale_shape_manual(values = c(21, 24)) +
  scale_colour_manual(values = c(plotColors[4], plotColors[3])) +
  scale_fill_manual(values = c(plotColors[1], plotColors[2])) +
  guides(color=guide_legend(override.aes=list(shape=21))) +
  guides(fill=guide_legend(override.aes=list(shape=21,fill=c(plotColors[1], plotColors[2]),color=c("white","white")))) +
  theme_minimal()

# save the PCA plot
ggsave("sample_pca.png", plot = sample_pca, device = "png")

# transpose of the transformed count matrix to get sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

# convert the distances to a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# update the column names
colnames(sampleDistMatrix) <- NULL

# create a list of continuous colors
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# view the similarities and dissimilarities between samples
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# store the clustering plot as a ggplot object
sample_clust <- as.ggplot(pheatmap(sampleDistMatrix,
                                   clustering_distance_rows=sampleDists,
                                   clustering_distance_cols=sampleDists,
                                   col=colors))

# save the plot to a png file
ggsave("sample_clustering.png", plot = sample_clust, bg = "white", device = "png")

##
# DE Analysis Contrasts
##

# standard differential expression analysis steps are wrapped into a single function
dds <- DESeq(dds)

# extract a results table with log2 fold changes, p values and adjusted p values
res <- results(dds)

# view the group names
resultsNames(dds)

# the effect of treatment in MANF (the main effect)
res_treatment_MANF <- results(dds, name="treatment_Food_Restriction_vs_Ad_lib")
#res_treatment_MANF <- results(dds, contrast=c("treatment","Food_Restriction","Ad_lib"))
# the effect of treatment in SARA
res_treatment_SARA <- results(dds, contrast=list(c("treatment_Food_Restriction_vs_Ad_lib", "treatmentFood_Restriction.mouselineSARA")))
# test if the treatment effect is different in SARA compared to MANF (interaction term)
res_interaction <- results(dds, name="treatmentFood_Restriction.mouselineSARA")
# the effect of mouseline in Ad_lib
res_mouseline <- results(dds, name="mouseline_SARA_vs_MANF")
# what is the difference between SARA and MANF without treatment?
#res_mouseline <- results(dds, contrast=c("mouseline", "SARA", "MANF"))
# with treatment, what is the difference between SARA and MANF? (double check)
res_mouseline_treatment <- results(dds, contrast=list(c("mouseline_SARA_vs_MANF", "treatment_Food_Restriction_vs_Ad_lib")))

# summarize some basic tallies
summary(res_treatment_MANF)
summary(res_treatment_SARA)
summary(res_interaction)
summary(res_mouseline)
summary(res_mouseline_treatment)

# retrieve information about which variables and tests were used
#mcols(res_mouseline_treatment_diff)$description

# order our results table by the smallest p value
resOrdered_treatment_MANF <- as.data.frame(res_treatment_MANF[order(res_treatment_MANF$pvalue),])
resOrdered_treatment_SARA <- as.data.frame(res_treatment_SARA[order(res_treatment_SARA$pvalue),])
resOrdered_interaction <- as.data.frame(res_interaction[order(res_interaction$pvalue),])
resOrdered_mouseline <- as.data.frame(res_mouseline[order(res_mouseline$pvalue),])
resOrdered_mouseline_treatment <- as.data.frame(res_mouseline_treatment[order(res_mouseline_treatment$pvalue),])

# add geneID column
resOrdered_treatment_MANF$geneID <- rownames(resOrdered_treatment_MANF)
resOrdered_treatment_SARA$geneID <- rownames(resOrdered_treatment_SARA)
resOrdered_interaction$geneID <- rownames(resOrdered_interaction)
resOrdered_mouseline$geneID <- rownames(resOrdered_mouseline)
resOrdered_mouseline_treatment$geneID <- rownames(resOrdered_mouseline_treatment)

# re-order columns
resOrdered_treatment_MANF <- resOrdered_treatment_MANF[, c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
resOrdered_treatment_SARA <- resOrdered_treatment_SARA[, c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
resOrdered_interaction <- resOrdered_interaction[, c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
resOrdered_mouseline <- resOrdered_mouseline[, c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
resOrdered_mouseline_treatment <- resOrdered_mouseline_treatment[, c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# save the ordered results to a csv file
write.csv(resOrdered_treatment_MANF, file="resOrdered_treatment_MANF.csv", quote = FALSE, row.names = FALSE)
write.csv(resOrdered_treatment_SARA, file="resOrdered_treatment_SARA.csv", quote = FALSE, row.names = FALSE)
write.csv(resOrdered_interaction, file="resOrdered_interaction.csv", quote = FALSE, row.names = FALSE)
write.csv(resOrdered_mouseline, file="resOrdered_mouseline.csv", quote = FALSE, row.names = FALSE)
write.csv(resOrdered_mouseline_treatment, file="resOrdered_mouseline_treatment.csv", quote = FALSE, row.names = FALSE)

# check the number of adjusted p-values were less than 0.05
sum(resOrdered_treatment_MANF$padj < 0.05, na.rm=TRUE)
sum(resOrdered_treatment_SARA$padj < 0.05, na.rm=TRUE)
sum(resOrdered_interaction$padj < 0.05, na.rm=TRUE)
sum(resOrdered_mouseline$padj < 0.05, na.rm=TRUE)
sum(resOrdered_mouseline_treatment$padj < 0.05, na.rm=TRUE)

# filter results by p-value
res05_treatment_MANF <- resOrdered_treatment_MANF[resOrdered_treatment_MANF$padj < 0.05,]
res05_treatment_SARA <- resOrdered_treatment_MANF[resOrdered_treatment_SARA$padj < 0.05,]
res05_interaction <- resOrdered_treatment_MANF[resOrdered_interaction$padj < 0.05,]
res05_mouseline <- resOrdered_treatment_MANF[resOrdered_mouseline$padj < 0.05,]
res05_mouseline_treatment <- resOrdered_treatment_MANF[resOrdered_mouseline_treatment$padj < 0.05,]

# remove NAs
res05_treatment_MANF <- na.omit(res05_treatment_MANF)
res05_treatment_SARA <- na.omit(res05_treatment_SARA)
res05_interaction <- na.omit(res05_interaction)
res05_mouseline <- na.omit(res05_mouseline)
res05_mouseline_treatment <- na.omit(res05_mouseline_treatment)

# test with the adjusted p-value cut offs
#res05_treatment_MANF <- results(dds, name="treatment_Food_Restriction_vs_Ad_lib", alpha=cutFDR, lfcThreshold=cutLFC)
#res05_treatment_SARA <- results(dds, contrast=list(c("treatment_Food_Restriction_vs_Ad_lib", "treatmentFood_Restriction.mouselineSARA")), alpha=cutFDR, lfcThreshold=cutLFC)
#res05_interaction <- results(dds, name="treatmentFood_Restriction.mouselineSARA", alpha=cutFDR, lfcThreshold=cutLFC)
#res05_mouseline <- results(dds, name="mouseline_SARA_vs_MANF", alpha=cutFDR, lfcThreshold=cutLFC)
#res05_mouseline_treatment <- results(dds, contrast=list(c("mouseline_SARA_vs_MANF", "treatment_Food_Restriction_vs_Ad_lib")), alpha=cutFDR, lfcThreshold=cutLFC)

# summarize the results
#summary(res05_treatment_MANF)
#summary(res05_treatment_SARA)
#summary(res05_interaction)
#summary(res05_mouseline)
#summary(res05_mouseline_treatment)

# save the filtered results to a csv file
write.csv(res05_treatment_MANF, file="treatment_MANF_padj0.05_lfc0.csv", quote = FALSE, row.names = FALSE)
write.csv(res05_treatment_SARA, file="treatment_SARA_padj0.05_lfc0.csv", quote = FALSE, row.names = FALSE)
write.csv(res05_interaction, file="interaction_padj0.05_lfc0.csv", quote = FALSE, row.names = FALSE)
write.csv(res05_mouseline, file="mouseline_padj0.05_lfc0.csv", quote = FALSE, row.names = FALSE)
write.csv(res05_mouseline_treatment, file="mouseline_treatment_padj0.05_lfc0.csv", quote = FALSE, row.names = FALSE)

##
# Results Exploration
##

# extract vst counts
vsdCounts <- assay(vsd)

# subset the vst data by the DGE set
vsdSubset_treatment_MANF <- vsdCounts[rownames(vsdCounts) %in% rownames(res05_treatment_MANF), ]
vsdSubset_treatment_SARA <- vsdCounts[rownames(vsdCounts) %in% rownames(res05_treatment_SARA), ]
vsdSubset_interaction <- vsdCounts[rownames(vsdCounts) %in% rownames(res05_interaction), ]
vsdSubset_mouseline <- vsdCounts[rownames(vsdCounts) %in% rownames(res05_mouseline), ]
vsdSubset_mouseline_treatment <- vsdCounts[rownames(vsdCounts) %in% rownames(res05_mouseline_treatment), ]

# setup annotation names
ann_groups <- data.frame(
  treatment = as.factor(factors$treatment),
  mouseline = as.factor(factors$mouseline),
  sex = as.factor(factors$sex),
  row.names = colnames(vsd))

# setup annotation colors
ann_colors <- list(
  treatment = c(plotColors[3], plotColors[4]),
  mouseline = c("black", "grey"),
  sex = c(plotColors[11], plotColors[6])
)

# setup annotation color names
names(ann_colors[["treatment"]]) <- unique(factors$treatment)
names(ann_colors[["mouseline"]]) <- unique(factors$mouseline)
names(ann_colors[["sex"]]) <- unique(factors$sex)

# create heatmap for DGE
treatment_MANF_vst_dge <- as.ggplot(
  pheatmap(vsdSubset_treatment_MANF, 
           scale="row", 
           annotation_col = ann_groups, 
           annotation_colors = ann_colors,
           #main="Heatmap of DE Genes", 
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12,
           color = colorRampPalette(c(plotColors[5], "white", plotColors[2]))(100))
)
# save the plot to a png file
ggsave("treatment_MANF_vst_dge.png", plot = treatment_MANF_vst_dge, bg = "white", device = "png", width = 9, height = 8, units = "in")

# create heatmap for DGE
treatment_SARA_vst_dge <- as.ggplot(
  pheatmap(vsdSubset_treatment_SARA, 
           scale="row", 
           annotation_col = ann_groups, 
           annotation_colors = ann_colors,
           #main="Heatmap of DE Genes", 
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12,
           color = colorRampPalette(c(plotColors[5], "white", plotColors[2]))(100))
)
# save the plot to a png file
ggsave("treatment_SARA_vst_dge.png", plot = treatment_SARA_vst_dge, bg = "white", device = "png", width = 9, height = 8, units = "in")

# create heatmap for DGE
interaction_vst_dge <- as.ggplot(
  pheatmap(vsdSubset_interaction, 
           scale="row", 
           annotation_col = ann_groups, 
           annotation_colors = ann_colors,
           #main="Heatmap of DE Genes", 
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12,
           color = colorRampPalette(c(plotColors[5], "white", plotColors[2]))(100))
)
# save the plot to a png file
ggsave("interaction_vst_dge.png", plot = interaction_vst_dge, bg = "white", device = "png", width = 9, height = 8, units = "in")

# create heatmap for DGE
mouseline_vst_dge <- as.ggplot(
  pheatmap(vsdSubset_mouseline, 
           scale="row", 
           annotation_col = ann_groups, 
           annotation_colors = ann_colors,
           #main="Heatmap of DE Genes", 
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12,
           color = colorRampPalette(c(plotColors[5], "white", plotColors[2]))(100))
)
# save the plot to a png file
ggsave("mouseline_vst_dge.png", plot = mouseline_vst_dge, bg = "white", device = "png", width = 9, height = 8, units = "in")

# create heatmap for DGE
mouseline_treatment_vst_dge <- as.ggplot(
  pheatmap(vsdSubset_mouseline_treatment, 
           scale="row", 
           annotation_col = ann_groups, 
           annotation_colors = ann_colors,
           #main="Heatmap of DE Genes", 
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize = 12,
           color = colorRampPalette(c(plotColors[5], "white", plotColors[2]))(100))
)
# save the plot to a png file
ggsave("mouseline_treatment_vst_dge.png", plot = mouseline_treatment_vst_dge, bg = "white", device = "png", width = 9, height = 8, units = "in")

# sort results by gene name
treatment_SARA_dge_results_tbl <- resOrdered_treatment_SARA[order(resOrdered_treatment_SARA$geneID), ]
treatment_MANF_dge_results_tbl <- resOrdered_treatment_MANF[order(resOrdered_treatment_MANF$geneID), ]
interaction_dge_results_tbl <- resOrdered_interaction[order(resOrdered_interaction$geneID), ]
mouseline_dge_results_tbl <- resOrdered_mouseline[order(resOrdered_mouseline$geneID), ]
mouseline_treatment_dge_results_tbl <- resOrdered_mouseline_treatment[order(resOrdered_mouseline_treatment$geneID), ]

# store LFCs
gene_IDs <- treatment_SARA_dge_results_tbl$geneID
init_DE <- rep("None", nrow(treatment_SARA_dge_results_tbl))
DGE_log2FoldChange <- data.frame(
  geneID = gene_IDs,
  SARA_LFC = treatment_SARA_dge_results_tbl$log2FoldChange,
  MANF_LFC = treatment_MANF_dge_results_tbl$log2FoldChange,
  SARA_DE = init_DE,
  MANF_DE = init_DE,
  Response_SARA_MANF = init_DE
)

# indicate DE gene set
treatment_SARA_dge_sig_tbl.keep <- DGE_log2FoldChange$geneID %in% res05_treatment_SARA$geneID
DGE_log2FoldChange$SARA_DE[treatment_SARA_dge_sig_tbl.keep] <- "SARA"
treatment_MANF_dge_sig_tbl.keep <- DGE_log2FoldChange$geneID %in% res05_treatment_MANF$geneID
DGE_log2FoldChange$MANF_DE[treatment_MANF_dge_sig_tbl.keep] <- "MANF"

# treatment
# SARA_MANF
DGE_log2FoldChange$Response_SARA_MANF <- paste(DGE_log2FoldChange$SARA_DE, DGE_log2FoldChange$MANF_DE, sep = "_")
# clean up
DGE_log2FoldChange$Response_SARA_MANF <- gsub("None_None", "None", DGE_log2FoldChange$Response_SARA_MANF)
DGE_log2FoldChange$Response_SARA_MANF <- gsub("SARA_None", "SARA", DGE_log2FoldChange$Response_SARA_MANF)
DGE_log2FoldChange$Response_SARA_MANF <- gsub("None_MANF", "MANF", DGE_log2FoldChange$Response_SARA_MANF)
DGE_log2FoldChange$Response_SARA_MANF <- gsub("SARA_MANF", "Both", DGE_log2FoldChange$Response_SARA_MANF)

# turn on scientific notation
options(scipen = 0)

# SARA_MANF plots
# order data frames
#custom_order <- c("None", "SARA", "MANF", "Both")
custom_order <- unique(DGE_log2FoldChange$Response_SARA_MANF)
DGE_log2FoldChange_sorted <- DGE_log2FoldChange %>%
  arrange(match(Response_SARA_MANF, custom_order))
# scatter plots of log2 FC with correlations
jpeg("treatment_SARA_MANF_dge_LFC_cor.jpg")
ggplot(DGE_log2FoldChange_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response_SARA_MANF)) + 
  geom_point(aes(alpha = Response_SARA_MANF)) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  stat_cor(method = "pearson") +
  #scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75, "Both" = 0.75)) + 
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "Both" = 0.75)) + 
  #scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6], "Both" = plotColors[7])) +
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[6], "Both" = plotColors[7])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()
