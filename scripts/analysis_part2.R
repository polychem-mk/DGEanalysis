# Differential gene expression (DGE) analysis of RNA sequencing data

# Part 2.
#       Libraries
#       Load data and variables
#       Exploratory data analysis 
#         Counts data
#         Metadata
#         DESeqDataSet object 
#         Count normalization
#         Surrogate variable analysis (SVA)
#         Variance stabilization
#         Correlation Heat map 
#         Principal component analysis (PCA)

#             Libraries                                             -----------

library(dplyr)       # Data manipulation 
library(stringr)
library(ggplot2)     # Data visualization 
library(pheatmap)    # heatmap plots
library(vsn)         # to make meanSdPlot()

library(gt)          # to make gt_table 
library(patchwork)   # combine multiple plots into a single layout

# Bioconductor packages 
library(sva)           # surrogate variable analysis 
library(DESeq2)        # DE analysis

# The package "hexbin" is required for vsn meanSdPlot() 
# install.packages("hexbin")    

# The package "Hmisc" is required for ggplot stat_summary() 
# install.packages("Hmisc")     

# For parts 2 and 3 install/update Bioconductor packages if missing:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "vsn", "sva", "apeglm", "AnnotationDbi", "org.Hs.eg.db"))

#             Load data and variables                               -----------

# Data that was processed in Part 1: 
#  'metadata_processed' has the following factor variables:
#           "cancer_type":  HC, breast, hepatobiliary, CRC, GBM, lung, pancreas  
#                 "batch":  1 2 3 4 5 6   
#   "mutational_subclass":  EGFR, EGFR_MET, HER2, HER2_PIK3CA, KRAS, KRAS_MET, 
#                           MET, PIK3CA, Triple_Negative, wt
metadata_processed <- readRDS("data/metadata_processed.rds")

# 'counts_matrix' has only numeric data; column names are the same as row 
# names in the 'metadata_processed'; row names are ENS gene identifiers
counts_matrix <- readRDS("data/counts_matrix.rds")

 # A custom function that creates a heatmap and adds annotations:
source("scripts/functions/plot_heatmap.R")  

# Colors for plots:
plot_colors  = list(
  batch = c("1" = "#6c8db1",  "2" =  "#b16c8a" ,"3" =  "#3CB371", 
            "4" ="#87CEEB"  , "5" = "#FF8674" , "6" = "#ADFF2F"),
  
  mutational_subclass = c(HER2 = "#f7dc6f", wt = "#85929e",
                          Triple_Negative = "#6ca05d", HER2_PIK3CA = "#b7950b",
                          PIK3CA ="#db804f", KRAS = "#f941d7", MET = "#fbb7fb",
                          EGFR_MET ="#6c408c", EGFR =  "#c0f1ea",
                          KRAS_MET =   "#c7b9ff" ) ,
  
  cancer_type = c(breast = "#C11B17", hepatobiliary = "#5CB3FF",
                  CRC = "#a629f5", GBM = "#a7d170", lung = "#f78909", 
                  pancreas = "#0711fe", HC = "#474747"))


# Labels for plots:
labels_ct = list(HC = "HC (control)",
                 breast = "breast", 
                 CRC = "CRC", 
                 GBM = "GBM",
                 hepatobiliary = "hepatobiliary",
                 lung = "lung",
                 pancreas ="pancreas" )


#             Exploratory data analysis                             -----------

##                 Counts data                                      -----------

# 'counts_matrix' contains RNA-seq count data for 57,736 genes and 283 samples. 
dim(counts_matrix)

# The raw values in 'counts_matrix' represent the counts of RNA sequencing reads 
# corresponding to specific genes and proportional to gene expression in the
# sample;
counts_matrix[1:5, 1:5]  # first 5 rows for first 5 samples

# About 50% of these genes have 0 counts . 
rs <- rowSums(counts_matrix)  # total reads counts for each gene
mean(rs == 0)

# The distribution of gene counts for the first sample shows that a large
# proportion of genes have zero or low read counts. For better visualization,
# we also plot the distribution of non-zero counts on a logarithmic scale: a 
# small proportion of genes have a very large number of reads (from thousands 
# to hundreds of thousands).
p1 <- counts_matrix %>%    # Histogram of the read counts for 1 sample
  dplyr::select(X3.Breast.Her2.ampl ) %>%
  ggplot(aes(x = X3.Breast.Her2.ampl  )) +   
  geom_histogram(stat = "bin", 
                 bins = 100,
                 fill = "#1f618d") +
  xlab("Expression counts") +
  ylab("Number of genes") +     
  ggtitle("all counts") +
  theme_light() 

# Histogram of the read counts, greater than 0, for sample 1; the number
# of reads is on a logarithmic scale
p2 <- counts_matrix %>%
  dplyr::select(X3.Breast.Her2.ampl ) %>%
  filter(X3.Breast.Her2.ampl >0) %>%
  ggplot(aes(x = X3.Breast.Her2.ampl  )) +   
  geom_histogram(stat = "bin", 
                 bins = 100,
                 fill = "#1f618d") +
  xlab("Expression counts in log scale") +
  ylab("Number of genes") +     
  scale_x_log10() +
  ggtitle("non zero counts") +
  theme_light() +
  theme(plot.margin = margin(0, 0, 0, 30))

p1 + p2 + plot_annotation(title = "Histogram of count data for the first sample")

##                 Metadata                                         -------------

# In this part we use the processed data, metadata_processed, that was saved in
# Part 1 and has the factor variables:
#   cancer_type (HC, breast, hepatobiliary, CRC, GBM, lung, pancreas);
#   batch (1 2 3 4 5 6);
#   mutational_subclass (EGFR, EGFR_MET, HER2, HER2_PIK3CA, KRAS, KRAS_MET, MET, 
#                        PIK3CA, Triple_Negative, wt)

# The row names are the same as column names of the counts_matrix.
metadata_processed[1:5, ]  # metadata for the first five samples

# The following bar charts show the distribution of values for the cancer_types 
# variable, colored by mutation type and batch.

# These data contain 6 cancer types and control samples ("HC", from healthy 
# individuals). Each group contains more than 30 samples, except for 
# hepatobiliary type cancer.

# Bar plots of the  cancer types, colored by mutational subclass and batch:
p0 <- metadata_processed %>%
  ggplot() + 
  theme_light() +
  theme(axis.text.x = element_text( size = 8),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8)       )

p1 <- p0 +
  geom_bar(aes(y = cancer_type, fill = mutational_subclass), width = 0.65) +
  scale_fill_manual(values = plot_colors$mutational_subclass) +
  ggtitle("Mutational subclass") +
  scale_y_discrete(labels = labels_ct) +
  theme(axis.text.y = element_text(size = 12),
        legend.key.size = unit(0.3, "cm"))

p2 <- p0 +
  geom_bar(aes(y = cancer_type, fill = batch), width = 0.65) +
  scale_fill_manual(values = plot_colors$batch) +
  ggtitle("Batch") +
  theme(axis.text.y = element_blank())

p1 + p2 

# A large proportion of samples are wild type, which is found in all cancer_type  
# groups. The most common mutation in the data is KRAS, which is found in 
# "pancreas", "lung", "hepatobiliary" and "CRC" samples. Mutations such as 
# triple negative, PIK3CA, HER2_PIK3CA and HER2 are specific to the "breast", 
# and MET and KRAS_MET are only in the "lung" cancer type.

# The graph on the left shows that there are 6 batches: batch 3 and 4 cover all    
# cancer type categories; batch 2 covers all cancer type categories except breast 
# cancer; batch 1 is only GBM, batch 5 and 6 are missing from the control group.


##                 DESeqDataSet object                              ----------

# There are two main goals in this part of analysis: first, we need to find  
# sources of variability in the data to account for it when performing DGE 
# analysis; and second, we want to check if all samples are similar (correlated) 
# within cancer_type groups to detect potential outliers.

# To explore patterns in this data, we will use unsupervised learning methods. 

# The DESeqDataSetFromMatrix() function creates a DESeqDataSet object and brings 
# together all the components needed for further analysis: the raw count matrix, 
# metadata, and the design formula.

# It is unlikely that there is no batch effect. But for now,
# we keep all the samples and set the design formula to 'cancer_type' only.

# The column names of the count matrix must be the same and in the same order
# as the row names of the metadata table. 
all(rownames(metadata_processed) == colnames(counts_matrix))

 # creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = metadata_processed,
                              design = ~ cancer_type)
dds

##                 Count normalization                             ----------

# The next step is to normalize the RNA-seq counts. This will make the data for
# different samples comparable despite different sequencing depths [6, 7]. 
# The estimateSizeFactors() function adds size factors to the DESeq2 object,  
# which gives an estimate of the sequencing depth; the sizeFactors() function 
# then allows access to the size factors.
dds <- estimateSizeFactors(dds)    # add size factors to  the DESeq2 object

# Data frame with size factors for each sample and total counts (column sums of 
# the count data);  batch and sample name are added from the metadata.
size_factors <- data.frame(size_factors = sizeFactors(dds),
                           col_sums = colSums(counts(dds)),
                           batch = metadata_processed$batch,
                           sample = rownames(metadata_processed))

# The following two graphs show how the size factors vary across different
# batches.

# Scatterplot of size factors versus column sums
p1 <- ggplot(data = size_factors, 
             aes(x = size_factors, y = col_sums, colour = batch)) + 
  geom_point(size = 1.5) +
  geom_smooth(aes(x = size_factors, y = col_sums, colour = "lm"), 
              method = lm, se = FALSE, linewidth = 0.5) +
  labs(x = "sizeFactors", y = "colSums", 
       title = "Size factors vs. column sums colored by batch") +
  scale_color_manual(values = plot_colors$batch) +
  theme_light() +
  theme(legend.title = element_blank())

# Boxplots of size factors stratified by batch variable
p2 <- ggplot(data = size_factors, aes(x = batch, y =  size_factors)) + 
  geom_boxplot(fill = "#a9cce3") +
  labs(y = "sizeFactors", title = "Size factors by batch") +
  theme_light()

p1 + p2

# For this dataset, size factors range from 0.2 to 2.5, and one sample has
# a size factor of about 4, possibly due to the 6 batches and large number 
# of samples.
range(size_factors$size_factors[ -which.max(size_factors$size_factors) ])
size_factors[which.max(size_factors$size_factors), ]

##                 Surrogate variable analysis (SVA)                -----------

# To confirm the presence of the batch effect, we perform a surrogate variable  
# analysis. 

# SVA can be used to find hidden patterns in RNA-seq data and estimate surrogate
# variables. These surrogate variables can then be added to the model [8].
# However, in this project, we will only use SVA to demonstrate the presence of
# a batch effect.

# We run SVA on the subset of the count matrix for genes with row means 
# greater than 1. The proportion of these genes is about 19%.
normalized_counts <- counts(dds, normalized = TRUE) # counts scaled by size factors
ind <- rowMeans(normalized_counts) > 1    # indices for genes with row means > 1
mean(ind)

# Filter the matrix with normalized counts and keep genes with row average > 1
normalized_counts <- normalized_counts[ind, ]

# For SVA we build two models:
  # model that includes condition of interest
model <- model.matrix(~cancer_type, colData(dds)) 

 # model removes the condition of interest
model0 <- model.matrix(~1, colData(dds)) 

 # estimate surrogate variables
svseq <- svaseq(normalized_counts, model, model0, n.sv = 2)

# The estimated surrogate variables are in the "sv" output of svaseq(), and can 
# be plotted and colored by batch and mutation subclass.
head(svseq$sv)

 # data for SVA plots:
sv = data.frame(sv1 = svseq$sv[ , 1], 
                sv2 = svseq$sv[ , 2], 
                batch = dds$batch,
                mutational_subclass = dds$mutational_subclass)

p0 <- ggplot(data = sv, aes(sv1, sv2, colour = batch)) +
  theme_classic() +
  theme( legend.title = element_blank())

# SVA plot, colored by 'batch'
p1 <- p0 +    # SVA plot, colored by 'batch'
  geom_point(aes(colour = batch)) +
  scale_color_manual(values = plot_colors$batch) +
  ggtitle("Batch")

# SVA plot, colored by 'mutational_subclass'
p2 <- p0 +   # SVA plot, colored by 'mutational_subclass'
  geom_point(aes(colour = mutational_subclass)) +
  scale_color_manual(values = plot_colors$mutational_subclass)  +
  ggtitle("Mutational subclass")

p1 + p2 + plot_annotation(title = "Surrogate variable analysis plots")

# The graph on the left shows that batches 5 and 6 are different from the
# batches present in the control samples; batches 2, 3, and 4 overlap, but are
# still fairly distinct groups; and there are only two samples for batch 1, one
# of which is close to batch 3 and the other to batch 2.

# The graph on the right shows that SVA did not detect any groups among the 
# mutational subclasses: all mutations are mixed. Below we will use PCA and 
# batch stratification to find out if subclasses have any effect.

# Given the presence of the batch effect, to conduct the DGE analysis and to be 
# able to compare the cancer groups with the control group, we leave only the 
# batches included in the "HC" (control) cancer group, namely 2, 3 and 4.
# We also remove the sample with size factor ~4

# indices for samples to be excluded from further analysis: 
samples_rm <- c(rownames(size_factors[which.max(size_factors$size_factors), ]),
                rownames(size_factors[which(size_factors$batch %in% c("1", "5", "6")), ]) )

samples_rm <- which(rownames(metadata_processed) %in% samples_rm) 

# Create a new count matrix and metadata 
counts_filtered <- counts_matrix[ , -samples_rm]
metadata_filtered <- metadata_processed[ -samples_rm, ]

# Drop unused levels  
metadata_filtered$batch <- droplevels(metadata_filtered$batch)

##                 Variance stabilization                           -----------
 
# The variance of gene expression (read counts) tends to increase with
# increasing mean expression level and may influence the results of PCA or 
# clustering analyses. To address this issue, we also apply variance-stabilizing
# transformations.

# The DESeq2 package has several functions that perform these types of 
# transformations [9]. 

# (1 and 2) varianceStabilizingTransformation() and vst() (a wrapper for the 
#        first function) are functions that calculate a variance stabilizing 
#        transformation;
# (3) rlog() is for Regularized Log Transformation.

# This functions
# - apply transformations that are useful when checking for outliers or as 
#   input for machine learning methods.
# - normalizes with respect to library size. 
# - does not require estimating the size factors.

# * rlog() is more robust when the size factors vary widely. 
# * However, if the number of samples is 50 or more, rlog() will take a long
#   time [from functions manuals]

# Therefore for this analysis we will use vst() as it is faster in computation   
# time for 240 samples. The sample with the maximum size factor was dropped as 
# well as batches 1, 5 and 6, this should narrow the range of size factors.

# The main goal of the next part is to find outliers in the samples. We want    
# to find samples that are less correlated (or less than a certain threshold) 
# within their group (cancer type). 

 # Create a new DESeq2 object using only samples for the batches 2,3 and 4;
 # also add 'batch' to the design matrix.
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = metadata_filtered,
                              design = ~ batch + cancer_type)

# We also set 'blind' argument to  FALSE;
# From function manual: If many of genes have large differences in counts due 
#                      to the experimental design, it is important to set
#                      blind=FALSE for downstream analysis.
vsd <- vst(dds, blind = FALSE) # run variance stabilizing transformation     
vsd_mat <- assay(vsd)  # and extract values

# Now the normalized data looks like this for the first 5 samples:
vsd_mat[1:5, 1:5]

# The scatter plots of the first two samples show that after applying the  
# variance stabilizing transformation, the variance of gene counts decreased 
# for genes with larger gene counts.

# On the Mean-Standard deviation plot the standard deviation for all 
# samples reaches a plateau as the row mean increases.

# Create mean vs. standard deviation plot for vst() output
p1 <- ggplot(data = counts_matrix[ , 1:2],
            aes(x = X3.Breast.Her2.ampl, y = X8.Breast.WT )) +
  geom_point(size = 0.5, shape = 21, color = "#1f618d") + 
  labs(x = "sample 1", y = "sample 2", title = "Before transformation") +
  theme_light()

p2 <- ggplot(data = vsd_mat[ , 1:2],
             aes(x = X3.Breast.Her2.ampl, y = X8.Breast.WT )) +
  geom_point(size = 0.5, shape = 21, color = "#1f618d") + 
  labs(x = "sample 1", y = "sample 2", title = "After VST") +
  theme_light()

 # Create Mean-SD plot 
p3 <- meanSdPlot(vsd_mat, ranks = FALSE, plot = FALSE)

p3 <- p3$gg + 
  ggtitle("Row standard deviations vs.
          row means") + 
  ylim(c(0, 5)) +
  theme_light()

p1 + p2 + p3

##                 Correlation Heat map                             -----------

# A correlation heat map can give us an idea of how similar (or different) 
# samples in a group are, and allows us to quickly spot outliers.
cor_mat <- cor(vsd_mat)  # correlation matrix 

# First we will look at the "HC" samples (control in "cancer_type"), all of    
# which are wild type (mutation subclass "wt").

# indices for control samples
ind_HC <-  which(metadata_filtered$cancer_type %in% c( "HC")) 
cor_mat_HC = cor_mat[ind_HC, ind_HC]  # subset control samples

# The heat map for control samples has column annotations with batch number.
cor_mat_HC %>%
  as.matrix() %>%
  pheatmap(annotation_col =  dplyr::select(metadata_filtered, batch),
           annotation_colors = list(batch = plot_colors$batch[c("2", "3","4")]),
           show_rownames = FALSE,
           show_colnames = FALSE,
           main = "Correlation heatmap for HC (control) samples")

# Above, the SVA plot showed that about half of the samples from batch 2 are    
# close to both batches 3 and 4, but there are many samples that form a group 
# apart from the other batches. 

# In the heat map, the samples from batch 2 also formed two groups: the bottom 
# group is similar to almost all samples, while the top one has a lower 
# correlation with many samples.

# Let's find possible outliers: samples that have a correlation of less than 
# 0.8 with 4 or more samples.
low_cor_ind = sapply(1:ncol(cor_mat_HC), function(j){
  sum(cor_mat_HC[ , j] < 0.8) 
})

# Indices for samples that are less correlated with other samples in the HC group
ind_rm_HC = which(rownames(metadata_filtered) %in%  
                    colnames(cor_mat_HC)[which(low_cor_ind >= 4)])

metadata_filtered[ind_rm_HC, ]

# Next we do the same for each cancer type: 
cancer_types = setdiff(as.character(unique(metadata_processed$cancer_type)), "HC")

# plot_heatmap() is a custom function (defined in file 'plot_heatmap.R') that 
# filters cancer type and creates heatmaps; it adds column annotations 
# ("mutational_subclass") and row annotations ("batch").

for(i in cancer_types){
  plot_heatmap(cm = cor_mat,
               md = metadata_filtered,
               level_to_compare = i)
}

# We filter out outliers in the same way as for "HC" samples.

 # samples with low correlations:
low_cor <- sapply(cancer_types, function(i){
  
  ind = which(metadata_filtered$cancer_type == i)   # indices for cancer_type i
  
  cor_mat_i = cor_mat[ind, ind]  # take a subset from the correlation matrix
  
  # number of samples with a correlation less than 0.8
  low_cor_ind = sapply(1:ncol(cor_mat_i), function(j){
    sum(cor_mat_i[ , j] < 0.8)
  })
  
  # get sample names for samples that have a correlation less than 0.8 with 4 
  # or more samples.
  colnames(cor_mat_i)[which(low_cor_ind >= 4)]
})

# Indices for samples to be excluded from further analysis 
ind_rm = which(rownames(metadata_filtered) %in%  c(unlist(low_cor))   )
ind_rm = c(ind_rm, ind_rm_HC)

##                 Principal component analysis (PCA)               ----------

# Finally, we can run the principal component analysis, and to make it simpler,  
# we use the plotPCA() function and get the data for only the first two 
# principal components.

# The vast majority of genes have zero or low counts, therefore we perform PCA 
# on a subset of genes with higher variance across the samples.

# Indices for genes with higher variance:
ind_var <- sort(matrixStats::rowVars(vsd_mat[ , -c(ind_rm)]), decreasing = TRUE)[1:500]
ind_var <- which(rownames(counts_filtered) %in% names(ind_var))

vsd_sub <- vsd[ind_var , -ind_rm]  # subset of the counts data after VST

# PC1 and PC2 data for the counts after variance stabilizing transformation
pca <- vsd_sub %>%  
  DESeq2::plotPCA(intgroup = c("cancer_type", "mutational_subclass", "batch"),
                       pcsToUse = 1:2,
                       returnData = TRUE)

# PCA plots (PC1 vs. PC2):
# Color by cancer_type
p1 <- ggplot(data =  pca, aes(x = PC1, y = PC2, colour = cancer_type)) +
  geom_point( size = 2) +
  scale_color_manual(values = plot_colors$cancer_type,
                     labels = labels_ct) +
  ggtitle("Cancer type") +
  theme_classic() + 
  theme(legend.position = "bottom", legend.title = element_blank())

# Color by batch
p2 <-ggplot(data =  pca, aes(x = PC1, y = PC2, colour = batch)) +
  geom_point( size = 2) +
  scale_color_manual(values = plot_colors$batch) +
  ggtitle("Batch") +
  theme_classic() + 
  theme(legend.position = "inside",
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(colour = "gray", linewidth = 0.3),
        legend.title = element_blank(),
        legend.direction = "vertical")

# Color by mutational_subclass
p3 <- ggplot(data =  pca, aes(x = PC1, y = PC2, colour = mutational_subclass))  +
  geom_point(size = 2) +
  scale_color_manual(values = plot_colors$mutational_subclass) +
  ggtitle("Mutational subclass") +
  theme_classic() +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 5))

p1 + p2 + p3 + plot_annotation(title = "The first two principal components")

# As above (SVA and correlation heatmaps) this data groups more by batch than 
# by mutational_subclass 

# Distribution of mutation subclass within batches: batch 3 covers the most of 
# the mutational_subclass (except for MET and KRAS_MET)
ggplot(data = metadata_filtered, aes(x = mutational_subclass, fill = batch)) +
  geom_bar(position = "fill", width = 0.5) + 
  scale_fill_manual(values = plot_colors$batch) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

# The following lapply() function calculates PC1 and PC2 using the plotPCA() 
# function for each batch (2 to 4).

pca_batch <- lapply(c("2", "3", "4"), function(i){
  
  ind_batch = which(vsd_sub$batch  == i)  # indices for batch i
  
  suppressMessages({
    vsd_sub[, ind_batch] %>%    # PC1 and PC2 data 
      DESeq2::plotPCA(intgroup = c("cancer_type", "mutational_subclass"),
                      pcsToUse = 1:2,
                      returnData = TRUE) %>% 
      mutate(batch = i)
  })
  
}) %>% bind_rows()

# Although some mutations form groups, these groups overlap with each other,
# and wt is mixed with all mutations.

# PCA plots faceted by batch and colored by mutational_subclass
ggplot(data =  pca_batch, aes(x = PC1, y = PC2, colour = mutational_subclass))+
  geom_point( size = 3) +
  scale_color_manual(values = plot_colors$mutational_subclass) +
  ggtitle("PCA plot, by batch and mutational subclass") +
  facet_grid(. ~ batch, 
             labeller = as_labeller(c("2" = "batch 2", "3" = "batch 3", "4" = "batch 4")) ) +
  theme_light() +
  theme(legend.text = element_text(size = 7),
        legend.position = "bottom",
        strip.background =element_rect(fill="gray95"),
        strip.text = element_text(colour = 'gray35')) 
  
# PCA analysis can separate the control samples (HC) from the cancer samples,
# but in batches 3 and 4, the GBM cancer group is close to or even mixed with 
# the control samples, and in batch 2, the lung cancer group is slightly mixed
# with the control samples. Based on the PCA plots, similar gene profiles are 
# likely to be observed between pancreatic cancer (pancreas), colorectal 
# cancer (CRC), and breast cancer.

# PCA plots faceted by batch and colored by cancer type
ggplot(data =  pca_batch, aes(x = PC1, y = PC2, colour = cancer_type))+
  geom_point( size = 3) +
  scale_color_manual(values = plot_colors$cancer_type
                     , labels = labels_ct
                     ) +
  ggtitle("PCA plots, by batch and cancer type") +
  facet_grid(. ~ batch, 
             labeller = as_labeller(c("2" = "batch 2", "3" = "batch 3", "4" = "batch 4"))) +
  theme_light() +
  theme(legend.text = element_text(size = 7),
        legend.position = "bottom",
        strip.background =element_rect(fill="gray95"),
        strip.text = element_text(colour = 'gray35')) 

# save filtered metadata and count matrix
metadata_filtered <- metadata_filtered[-ind_rm, ]
counts_filtered <- counts_filtered[ , -ind_rm]

# saveRDS(metadata_filtered, "data/metadata_filtered.rds")
# saveRDS(counts_filtered, "data/counts_filtered.rds")



