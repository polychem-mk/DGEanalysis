# Part 3.
#       Libraries
#       Load data and variables
#       Differential expression analysis
#       Results
#         MA plot and Volcano plot
#         Annotation
#         plotCounts 
#         DEresults
#         pheatmap
#       Compare all results
#       Conclusion
#       References
#       sessionInfo

#                Libraries                                   ------------------

library(dplyr)  # Data manipulation 
library(tidyr)
library(stringr)

library(ggplot2)   # requires package Hmisc for stat_summary()
library(gt)        # to make gt_table 
library(patchwork) # combine multiple plots into a single layout

# Bioconductor packages 
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(pheatmap) # heatmap plots

# For parts 2 and 3 install/update Bioconductor packages if missing:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "vsn", "sva", "apeglm", "AnnotationDbi", "org.Hs.eg.db"))

#                Load data and variables                     ------------------

# Processed data was filtered in Part 2 and saved as rds files:
# metadata and count matrix.

# Metadata with the following factor variables:
#           "cancer_type":  HC, breast, hepatobiliary, CRC, GBM, lung, pancreas  
#                 "batch":  2 3 4   
#   "mutational_subclass":  EGFR, EGFR_MET, HER2, HER2_PIK3CA, KRAS, KRAS_MET, 
#                           MET, PIK3CA, Triple_Negative, wt

# Load metadata and counts matrix:
metadata_filtered <- readRDS("data/metadata_filtered.rds")

counts_filtered <- readRDS("data/counts_filtered.rds")

 # A custom function summarizes the results for DGE analysis and returns:
 # - table with the top n differentially expressed genes; 
 # - summary and keywords; 
 # - ggplot for normalized counts.
source("scripts/functions/DEresults.R")

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

#                Differential gene expression analysis       -----------------

# GSE68086 is the large dataset, containing more than two groups to compare.  
# One way is to use a subset for each cancer type and perform pairwise analyses. 
# Another way is to use the DESeq2 package, which supports multi-factor designs,
# multi-level factors and experiments with many samples [9, 10].

#  Creating the DESeq2 object              

# For DGE analysis, we create a new DESeqDataSet using the count matrix 
# 'counts_filtered' and the metadata table 'metadata_filtered' that were saved
# as RDS files in part 2 'analysis_part2.R'.

# From PCA and SVA we know there is a batch effect, so we include it in our  
# design formula, with 'cancer_type' last since that is the variable we are
# interested in. We also include 'mutational_subclass' even though PCA showed
# that all mutations are mixed with "wt", but we see distinct  groups of KRAS 
# and EGFR_MET mutations.

 # check reference level of the cancer_type
levels(metadata_filtered$cancer_type)

# Check the order
all(rownames(metadata_filtered) == colnames(counts_filtered))
 
# Create DESeqDataSet 
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = metadata_filtered,
                              design = ~ batch + mutational_subclass + cancer_type)
dds

#   Run DESeq()  

# The DESeq() function performs differential expression analysis. DESeq2 also 
# includes shrinkage of log2 fold changes to improve the accuracy needed for 
# gene visualization and ranking. There are two ways:
# (1) run DESeq() and perform shrinkage afterwards using function lfcShrink()
# (2) another way is to set the 'betaPrior' argument to TRUE [9].

de <- DESeq(dds, betaPrior = TRUE )    # Run analysis, ~5 minutes

# The plot of dispersion estimates shows that dispersion decreases as the mean 
# increases, and the raw dispersions follow the fitted line fairly well, 
# indicating a reasonable fit of the data to the DESeq2 model.
plotDispEsts(de, cex = 0.5, legend = FALSE)

#                Results                                     --------------

# We first look at the results for one cancer type, CRC (contrast set to 'CRC' 
# and 'HC'), and set alpha to 0.01 so that results for adjusted p-values less 
# than 0.01 are colored in the MA and Volcano plots.
res_CRC <- results(de,  # extract results for CRC           
                   contrast = c("cancer_type", "CRC", "HC"),
                   alpha = 0.01)

# There are a lot more down-regulated genes (19%) than up-regulated (~2%)
# (percentage of all genes in the data set)
summary(res_CRC)

##                 MA plot and Volcano plot                --------------

# There are two types of graphs that visualize the log2 fold change in gene 
# expression in one group of samples (cancer) compared to another (control). 
# The MA plot is a scatterplot with the log2 fold changes on the y-axis and the 
# mean of the normalized counts on the x-axis. The Volcano plot displays the 
# log2 fold changes on the x-axis against statistical significance (p-values or 
# adjusted p-values) on y-axis.

# Data for MA plot 
p1 = DESeq2::plotMA(res_CRC, returnData = TRUE)  

# ggplot for the for MS plot
p1 = p1 %>%                    
  ggplot(aes(x = mean + 0.001, # add 0.001 to avoid infinite values
             y = lfc,
             color = isDE)) +
  geom_point(size = 0.8, alpha = 0.8) +
  # The DESeq2::plotMA() function displays the results in logarithmic scale,  
  # but for the ggplot we need to add a logarithmic x-axis scale.
  scale_x_log10() +
  scale_color_manual(values = c("gray", "blue")) +
  guides(color = "none") +
  labs(x = "mean of normalized counts",
       y = "log fold change",
       title = "MA plot") +
  theme_light()

res_CRC <- res_CRC %>%       # data for Volcano plot
  as.data.frame() %>%
  mutate(threshold = padj < 0.01)

p2 <- res_CRC %>%            # ggplot for Volcano plot
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj),
             color = threshold)) +
  geom_point(size = 0.8, alpha = 0.8) +
  labs(title = "Volcano plot", 
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  guides(color = "none") +
  scale_color_manual(values = c( "gray", "blue")) +
  theme_light()

# The MA and Volcano plots show this trend: there are more values with negative  
# log fold change values corresponding to p-values less than 0.01.

# If we check the relevant publications, we can find articles that report that 
# more of differentially expressed genes are down-regulated in cancer cell 
# samples; however, the difference is not as significant as in these data[11,12]
p1 + p2

##                 Annotation                              ----------------

# Adding the gene name and gene type can help to understand which types of  
# genes are expressed differently and compare the results with literature data.

# For this purpose, we use the org.Hs.eg.db package to add "SYMBOL", "GENENAME" 
# and "GENETYPE" to the results table with top 20 differently expressed genes,
# by "ENSEMBL" (Ensembl id).
# Hs.eg.db is an organism specific annotation package for Human genome.

# We take results with adjusted p-values less than 0.01 and order them by the 
# absolute value of the 'stat' (Wald statistic), which is the log2FoldChange 
# divided by lfcSE, to find results corresponding to larger log fold changes 
# and smaller standard errors.

# add variable ENSEMBL which will be used to match records from org.Hs.eg.db
res_CRC$ENSEMBL <- rownames(res_CRC) 

res_CRC <- res_CRC %>%       # results for top 20 differently expressed genes
  filter(padj < 0.01) %>%
  arrange(desc(abs(stat))) %>%
  slice_head(n = 20)

ensembl_ids = res_CRC$ENSEMBL # Ensembl IDs for top 20 genes

# Annotation data for Ensembl identifiers in res_CRC
ann_data <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = ensembl_ids, 
                                  # columns from org.Hs.eg.db
                                  columns = c("SYMBOL", "GENENAME", "GENETYPE" ),
                                  keytype = "ENSEMBL")
res_CRC <- res_CRC %>%    # add symbol, gene name and gene type to 'res_CRC'
  left_join(ann_data, by = "ENSEMBL")

# Most of these genes are protein coding, and there are many ribosomal proteins.  
sum(str_detect(res_CRC$GENETYPE, "protein-coding"), na.rm = TRUE)/length(res_CRC$GENETYPE)

sum(str_detect(res_CRC$GENENAME, "ribosomal"), na.rm = TRUE)/length(res_CRC$GENETYPE)

##                 plotCounts                              --------------

# The plotCounts() function (DESeq2) plots the normalized counts for a single
# gene. We can look at one gene that is differentially expressed in CRC, but
# make a plot for all cancer types.

ind <- res_CRC$ENSEMBL[2] # index for one gene from the results table for CRC

plot_title <- paste("Normalized counts of ", ind, " gene, by batch" , sep = "")

# Plot  normalized counts for one gene for all groups in the cancer_type 
# and color by batch
set.seed(1)
plotCounts(dds,
           ind, 
           intgroup = "cancer_type",
           returnData = TRUE)  %>%
  bind_cols(batch = metadata_filtered$batch,
            mutational_subclass = metadata_filtered$mutational_subclass) %>%
  ggplot(aes(x = cancer_type, y = count, colour = batch)) +
  geom_jitter(width = 0.3, size = 2, shape = 16, alpha = 0.5) +
  stat_summary(fun.data = mean_sdl,
               position = position_dodge(width = 0.5),
               linewidth = 1.1) +
  # The plotCounts() function displays the results in logarithmic scale,  
  # but for the ggplot we need to add a logarithmic y-axis scale.
  scale_y_log10() +
  scale_color_manual(values = plot_colors$batch) +
  scale_x_discrete(labels = labels_ct) +
  ggtitle(plot_title) + theme_light() + 
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank()) 

# Overall, batch 2 has lower means and higher standard deviations. Batches 3 and
# 4 are more similar, however, in the PCA plot, batch 4 showed better separation
# of cancer types, so for further analysis we will focus on the results for
# samples from batch 4.

# The following graph is a plot of normalized counts for the top 20
# differentially expressed genes in CRC group (Ensembl IDs have been
# replaced with gene symbols (SYMBOL)).

# Get data using plotCounts() function for all 20 genes in the results table
# for CRC:
res_CRC_plot <- lapply(res_CRC$ENSEMBL, function(i){
  
  plotCounts(dds, i, intgroup = "cancer_type",returnData = TRUE)  %>%
    # add batch and mutational_subclass variables from metadata_filtered:
    bind_cols(batch = metadata_filtered$batch,
              mutational_subclass = metadata_filtered$mutational_subclass) %>%
    filter(cancer_type %in% c( "CRC", "HC"), batch == "4") %>%
    mutate(ENSEMBL = i)
  
}) %>% bind_rows()

# Add Ensembl IDs:
res_CRC_plot <- res_CRC_plot %>%
  left_join(res_CRC[ , c("ENSEMBL", "SYMBOL")]) %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL), ENSEMBL , SYMBOL))

# Make a plot:
set.seed(1)
ggplot(res_CRC_plot, aes(x = SYMBOL , y = count, colour = cancer_type)) +
  geom_jitter(width = 0.3, size = 2, shape = 16, alpha = 0.5) +
  stat_summary(fun.data = mean_sdl,
               position = position_dodge(width = 0.2),
               linewidth = 1.1) +
  scale_y_log10() +
  scale_color_manual(values = c("#5d6d7e" , "#F08080" ),
                     labels = c("HC (control)", "CRC") ) +
  labs(x = "" , y = "normalized counts", title = "Cancer type CRC, top 20 genes") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())

# We can compare these results with the available publications, keeping in mind
# that the GSE68086 RNA-seq data are obtained from blood samples, whereas many 
# available publications typically study tumor tissue.

# For example, RP genes (ribosomal proteins) are significantly down-regulated
# in CRC cancer samples for this dataset. There are studies showing increased
# expression of RP genes in CRC samples[13], while other studies found that RP
# genes are down-regulated[14,15].

res_CRC %>%   # a subset of the results table that includes only RP genes
  dplyr::select(log2FoldChange, SYMBOL, GENENAME) %>%
  arrange(SYMBOL) %>%
  filter(str_detect(SYMBOL, "RP"))

# Even greater log-fold change in ANKR (ankyrin repeat domain) gene expression.
# And we can find articles reporting that ANKR (ankyrin repeat domain) mRNA
# expression is significantly lower in CRC tumor tissues than in control
# samples[16].
res_CRC %>%    # a subset of the results table that includes only ankyrin genes
  dplyr::select(log2FoldChange, SYMBOL, GENENAME) %>%
  arrange(SYMBOL) %>%
  filter(str_detect(SYMBOL, "ANKR"))

##                 DEresults                               -------------

# Now let’s look at other cancer types, but put the above code into a custom
# function DEresults(), that returns a list of 3:

#  1 Table with top n differentially expressed genes.
#  2 Summary and keywords for these genes. 
#  3 ggplot for normalized counts filtered for the specified batch (or all
#    batches if ‘batch=NULL’).

# The DEresults() function shows the results for all batches. It does not rerun
# the DESeq() function, but only retrieves the results. However, the data for
# the normalized counts plot can be filtered using the ‘batch’ argument.

# Get top 20 results for lung cancer type using custom function DEresults()
DEresults_lung <- DEresults(
  de = de,                  # output of DESeq() function
  md = metadata_filtered ,  # metadata
  level.1 = "HC",     # the base level (control) and  
  level.2 = "lung",   # the level to compare in the factor variable 'cancer_type'
  alpha = 0.01,       # cut off for padj
  batch = 4,          # batch number or NULL
  sort.by = "stat",   # sort results by "stat", "padj"  or "log2FoldChange"
  n = 20              # number of results
  )

# The lung cancer results also indicate that the vast majority of differentially
# expressed genes are down-regulated and protein coding.

DEresults_lung[[3]]   # normalized counts plot 
DEresults_lung[[2]]   # summary and keywords 

# About 41% of the genes in GSE68086 do not have a corresponding GENETYPE in the 
# org.Hs.eg.db database. But we can compare the distribution of the available gene
# types in the entire dataset and those genes that are differentially expressed.

# Data frame for all gene types corresponding to genes in the GSE68086 dataset
gene_types_all <- AnnotationDbi::select(org.Hs.eg.db, 
                      keys = rownames(counts_filtered), 
                      columns = "GENETYPE" ,
                      keytype = "ENSEMBL")  %>%
  # org.Hs.eg.db may have more than one row per Ensembl ID, so keep only
  # the first match.
  distinct(ENSEMBL, .keep_all = TRUE)

mean(is.na(gene_types_all$GENETYPE))  # proportion of the missing values

# Drop NAs in 'gene_types_all' data frame and count gene types:
gene_types_all <- gene_types_all %>%
  drop_na() %>%
  group_by(GENETYPE) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  mutate(prop = n/sum(n))

# The org.Hs.eg.db contains about 55.8% protein-coding genes (shown in gray),
# while differentially expressed genes are mainly protein-coding genes: 85% for 
# CRC samples and more than 90% for lung cancer samples (shown in blue).

gene_types_all

# "DEresults_lung[[2]][[2]]" is a table that contains the number of gene types 
# for genes with padj values less than 0.01 (alpha)
gene_types_de <- DEresults_lung[[2]][[2]] %>%  
  # n is the number of gene types regardless of up or down regulation:
  mutate(n = up + down) 

# Plot: proportions of gene types in org.Hs.eg.db and differentially expressed
# genes in lung cancer type
ggplot(data = gene_types_all) + 
  geom_col(aes(x = GENETYPE, y = n/sum(n)),
           fill = "gray", width = 0.4, just = 0.5) +
  geom_col(data = gene_types_de,
           aes(x = gene.type, y = n/sum(n)),
           fill = "steelblue", width = 0.4, just = -0.5) +
  labs(title = "Proportion of gene types in org.Hs.eg.db \n(gray) compared to proportion of \ndifferentially expressed genes (blue)\nin lung cancer samples",
       y = "proportion of genes",
       x = "gene type") +
  guides(fill = "none") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

##                 pheatmap                                ----------

## The heat map allows us to see the bigger picture. For better visualization,  
# we take samples from batch 4 and the genes with the most pronounced 
# differential expression across all cancer types.

# Levels of the variable cancer_type, except for "HC"
cancer_types = setdiff(as.character(unique(metadata_filtered$cancer_type)), "HC")

# Get top 30 results for all cancer_type groups
DEresults_all <- lapply(cancer_types, function(i){  # takes ~1 minute
  DEresults(de = de,
            md = metadata_filtered ,
            level.1 = "HC",
            level.2 = i,
            alpha = 0.01,
            batch = 4,
            n = 30  )
})

# Unique Ensembl IDs from 'DEresults_all':
gene_ids <- lapply(1:6, function(i){
  DEresults_all[[i]][[1]]$ENSEMBL
}) %>% unlist() %>% unique()

batch_ind <- which(metadata_filtered$batch == "4") # indices for batch 4

# note that there are only 103 unique genes, many genes are common across  
# several or all cancer types
length(gene_ids)

# We take a subset of counts after the variance stabilizing transformation that
# includes samples from batch 4 and the 103 differentially expressed genes 
# identified above. To create a more uniform heat map, we subtract the row means
# from the resulting matrix.

# Apply variance stabilizing transformation
vst <- varianceStabilizingTransformation(dds) %>% assay()

# Filter normalized count matrix for batch 4 and differentially expressed genes
vst <- vst[gene_ids, batch_ind] 
vst <- vst - rowMeans(vst)

# Data frame with metadata for the batch 4:
metadata_b4 = metadata_filtered[batch_ind, ]

# Make a heat map
pheatmap(vst, 
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = dplyr::select(metadata_b4, cancer_type),
         annotation_colors = plot_colors)

# In this heat map the samples are in the columns and the genes are in the rows

# The controls ("HC") are well defined, but the "GBM" samples are closer to them 
# and even mixed with them. Thus, it is likely more difficult to identify genes 
# that will differ significantly between the "GBM" and controls. In addition, 
# the cancer types are quite mixed with each other, although lung , CRC, and 
# breast cancer tend to cluster more.

# Normalized counts of top genes overlap more in case of GBM samples for batch 4
DEresults_all[[which(cancer_types == "GBM")]][[3]]
# (a better separation in normalized counts if we make the same graph for batch 2)

# Compared with CRC and lung cancer, the ratio of down to up regulated genes is
# not as large in GBM. More over, most of the top 20 genes are up-regulated.
DEresults_all[[which(cancer_types == "GBM")]][[2]][[1]]

#                Compare all results                        ---------

# The 'DEresults_all' list contains results for the six cancer types.
# The following summary shows the percentage of genes that are up- and 
# down-regulated. Given that we found a lot of RP and ANKRD genes, let's add how 
# many of these genes are present in the top 30 results for all cancer types.
# We can also look at the main keywords extracted from the gene names.
lapply(1:6, function(i){
  
  # The second element of the DEresults() output contains a summary of all
  # differentially expressed genes (not just the top results) that correspond to 
  # results with a p-value less than 0.01.
  up = sum(DEresults_all[[i]][[2]][[2]]$up)    # number of up-regulated genes 
  down = sum(DEresults_all[[i]][[2]][[2]]$down) # number of down-regulated genes 
  n = up+down
  
  # The first element of the DEresults() output is a table with the top n 
  # results, so we look at the number of RP and ANKRD genes in that table:
  ribosomal_protein = sum(str_detect(DEresults_all[[i]][[1]]$GENENAME,
                                     "ribosomal protein"),
                          na.rm = TRUE)
  ankyrin = sum(str_detect(DEresults_all[[i]][[1]]$GENENAME,
                           "ankyrin"),
                na.rm = TRUE)
  
  data.frame(cancer_type = cancer_types[i],
             up = round(100*up/n, 2),
             down = round(100*down/n, 2) ,
             ribosomal_protein = ribosomal_protein,
             ankyrin = ankyrin,
             keywords = toString(DEresults_all[[i]][[2]][[3]][1:10]))
  
}) %>% bind_rows()

# Most of the differentially expressed genes found in this analysis are 
# down-regulated.

# Some top differentially expressed genes are found in more than one cancer type.

 # Calculate the frequency of occurrence of the top genes in all types of cancer:
gene_by_ct <- sapply(1:6, function(i){
  
  # Ensembl IDs (ENSEMBL) for top genes are in 'gene_ids'
  sapply(gene_ids, function(j){
    as.numeric(j %in% DEresults_all[[i]][[1]]$ENSEMBL  )
    
  }) 
}) %>% as.data.frame()

# The row sums shows how many times each gene from the 'gene_ids' vector occurred 
# across all cancer types
rs <- rowSums(gene_by_ct)

# Add column names to 'gene_by_ct': order of cancer_type in 'DEresults_all' is
# same as in 'cancer_types' defined above
colnames(gene_by_ct) <- cancer_types
 
gene_by_ct$rs <- rs  # add row sums variable

gene_by_ct$ENSEMBL <- rownames(gene_by_ct) # add ENSEMBL variable

 # Annotation data for Ensembl identifiers in res_CRC
ann_data <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = gene_by_ct$ENSEMBL, 
                                  # columns from org.Hs.eg.db
                                  columns = c("SYMBOL", "GENENAME", "GENETYPE"),
                                  keytype = "ENSEMBL") %>%
  distinct(ENSEMBL, .keep_all = TRUE)  # keep the first match for Ensembl ID

# Add symbol, gene name and gene type to 'gene_by_ct'
gene_by_ct <- gene_by_ct %>%    
  left_join(ann_data, by = "ENSEMBL") %>%
  # if SYMBOL is missing, replace it with ANSEMBL
  mutate(SYMBOL = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

# The table 'gene_by_ct' shows that GBM has only one gene in common with other
# cancer types, while the other 5 cancer types have many genes in common. The
# ANKRD and RP genes are again present in this table.
gene_by_ct %>% 
  filter(rs >=4) %>%      # filter genes common to at least 4 cancer types
  dplyr::select(SYMBOL,GENENAME, breast, hepatobiliary,CRC,GBM,lung, pancreas) %>%
  arrange(SYMBOL)

# the following genes appeared only for one cancer type:
sapply(1:6, function(i){
  genes = gene_by_ct$SYMBOL[which(gene_by_ct$rs == 1 & gene_by_ct[, i] == 1)]
  paste(cancer_types[i], ": ", toString(sort(genes)), sep = "")
})

# We can visualize 'gene_by_ct' using a correlation heat map.

# Get numeric columns from 'gene_by_ct' and convert them to matrix
n_common_genes <- as.matrix(gene_by_ct[ , 1:6] == 1)  

# Determine the number of common genes across all cancer types (from top 30 
# results)
n_common_genes <- crossprod(n_common_genes) 
n_common_genes

# The following heat map shows the number of common gens and correlations 
# (as color) between all cancer types 
pheatmap(cor(gene_by_ct[,1:6]),     # create heat map
         display_numbers = n_common_genes )

# In terms of gene expression, these data suggest that the pancreas and CRC  
# samples are similar and GBM is the least similar to other cancer types in  
# this dataset

# We can set 'level.1' to another cancer type in the DEresults() function and 
# compare the normalized counts for two cancer types. 

# The following graph shows that there are very few differentially expressed
# genes for CRC and pancreas and the counts overlap, even when setting alpha
# to 0.05.
DEresults_CRC_panc <- DEresults(de = de,
                                md = metadata_filtered ,
                                level.1 = "CRC",   
                                level.2 = "pancreas",
                                alpha = 0.05,
                                sort.by = "log2FoldChange",
                                batch = 4,
                                n = 20 )
DEresults_CRC_panc[[3]]

# We can get better results for hepatobiliary and lung because there are much
# fewer common differentially expressed genes.
DEresults_hep_lung <- DEresults(de = de,
                                  md = metadata_filtered ,
                                  level.1 = "hepatobiliary",   
                                  level.2 = "lung",
                                  alpha = 0.01,
                                  sort.by = "log2FoldChange",
                                  batch = 4,
                                  n = 20 )

DEresults_hep_lung[[3]]

#                Conclusion                                 ------------

# The GSE68086 datasets contain counts of RNA reads for 285 samples and 57,736 
# genes. These data are from a high-throughput sequencing experiment and include 
# 55 samples from healthy individuals (HC, control group) and 228 samples from 
# patients with six different cancer types: lung cancer, colorectal cancer, 
# pancreatic cancer, breast cancer, glioblastoma, and hepatobiliary carcinoma.

# • The main variables that can be extracted from the raw metadata file are 
# batch number (batch), cancer type (cancer_type), and mutation subclass 
# (mutational_subclass).

# ◦ 61.1% of the samples are wild type (wt), which is found in all cancer groups. 
# The most common mutation in the data is the KRAS (22.3%), which is found the
# pancreas, lung, hepatobiliary, and CRC samples. Mutations such as 
# triple_negative, PIK3CA, HER2_PIK3CA and HER2 are specific to breast cancer,
# and MET and KRAS_MET are found only in the lung cancer group.

# ◦ All samples were obtained in 6 batches, however only batches 2, 3 and 4 are
# present in the control group (HC).

# ◦ During the cleaning/filtering process, 53 samples were discarded, including
# two samples with low read counts, samples that were less correlated within
# their cancer_type group, and batches 1, 5, and 6.

# • The counts_matrix that was used as input to the DESeq() function 
# contained the raw counts of RNA sequencing reads corresponding to 57,736 genes
# and 232 samples.

# ◦ About 50% of these genes have 0 raw counts

# ◦ To apply unsupervised machine learning methods, the counts were normalized
# and a variance stabilizing transformation was applied (DESeq2 vst() function  
# was used for the entire dataset and varianceStabilizingTransformation() for a
# subset).

# • Surrogate variable analysis (SVA) and Principal component analysis (PCA)
# showed the presence of the batch effect.

# • We observe an increase in the proportion of protein-coding genes among 
# differentially expressed genes.

# • For all cancer types, there are more down-regulated genes than up-regulated
# genes, but for GBM this difference is smaller than for the other 5 cancer
# types in this dataset.

# • If we take the 30 genes with the highest differential expression for each
# cancer type, we get a total of 103 unique genes, since there are genes that
# are common to at least two cancer types.
 
# • For CRC, lung and pancreas, about half of the genes in the top results are
# RP genes.

# • There are 24 out of 30 same genes in top results for CRC and pancreas.
 
# • Gene expression in GBM samples differed least from control samples among all
#   cancer types in this dataset.

#                References                                 -----------

# (1) Li, M.; Sun, Q.; Wang, X. Transcriptional Landscape of Human Cancers.
# Oncotarget 2017, 8 (21), 34534–34551. https://doi.org/10.18632/oncotarget.15837.
# (2) Xu, H.; Ma, Y.; Zhang, J.; Gu, J.; Jing, X.; Lu, S.; Fu, S.; Huo, J. 
# Identification and Verification of Core Genes in Colorectal Cancer. 
# BioMed Research International 2020, 2020 (1), 8082697.
# https://doi.org/10.1155/2020/8082697.
# (3) Hu, H.; He, J.; Zhao, H. Construction and Validation of a Diagnostic Model
# for Cholangiocarcinoma Based on Tumor-Educated Platelet RNA Expression Profiles.
# Oncologie 2025, 27 (2), 277–293. https://doi.org/10.1515/oncologie-2024-0520.
# (4) Ge, X.; Yuan, L.; Cheng, B.; Dai, K. Identification of Seven 
# Tumor-Educated Platelets RNAs for Cancer Diagnosis. Journal of Clinical
# Laboratory Analysis 2021, 35 (6), e23791. https://doi.org/10.1002/jcla.23791.
# (5) Best, M. G.; Sol, N.; Kooi, I.; Tannous, J.; Westerman, B. A.; Rustenburg,
# F.; Schellen, P.; Verschueren, H.; Post, E.; Koster, J.; Ylstra, B.; Ameziane, 
# N.; Dorsman, J.; Smit, E. F.; Verheul, H. M.; Noske, D. P.; Reijneveld, J. C.;
# Nilsson, R. J. A.; Tannous, B. A.; Wesseling, P.; Wurdinger, T. RNA-Seq of
# Tumor-Educated Platelets Enables Blood-Based Pan-Cancer, Multiclass, and 
# Molecular Pathway Cancer Diagnostics. Cancer Cell 2015, 28 (5), 666–676. 
# https://doi.org/10.1016/j.ccell.2015.09.018.
# (6) Wang, L.; Felts, S. J.; Van Keulen, V. P.; Pease, L. R.; Zhang, Y. 
# Exploring the Effect of Library Preparation on RNA Sequencing Experiments.
# Genomics 2019, 111 (6), 1752–1759. https://doi.org/10.1016/j.ygeno.2018.11.030.
# (7) Evans, C.; Hardin, J.; Stoebel, D. M. Selecting Between-Sample RNA-Seq
# Normalization Methods from the Perspective of Their Assumptions. Briefings in
# Bioinformatics 2017, 19 (5), 776–792. https://doi.org/10.1093/bib/bbx008.
# (8) Love, M. I.; Huber, W.; Anders, S. Moderated Estimation of Fold Change and 
# Dispersion for RNA-Seq Data with DESeq2. Genome Biology 2014, 15 (12), 550. 
# https://doi.org/10.1186/s13059-014-0550-8.
# (9) Love, M. I.; Anders, S.; Huber, W. Analyzing RNA-Seq Data with DESeq2, 2025.
# (10) Danielsson, F.; Skogs, M.; Huss, M.; Rexhepaj, E.; O’Hurley, G.;
# Klevebring, D.; Pontén, F.; Gad, A. K. B.; Uhlén, M.; Lundberg, E. Majority of
# Differentially Expressed Genes Are down-Regulated During Malignant 
# Transformation in a Four-Stage Model. Proceedings of the National Academy of 
# Sciences 2013, 110 (17), 6853–6858. https://doi.org/10.1073/pnas.1216436110.
# (11) Long, T.; Liu, Z.; Zhou, X.; Yu, S.; Tian, H.; Bao, Y. Identification of
# Differentially Expressed Genes and Enriched Pathways in Lung Cancer Using 
# Bioinformatics. Molecular Medicine Reports 2019, 19 (3), 2029–2040. 
# https://doi.org/10.3892/mmr.2019.9878.
# (12) Jing, X. Ribosomal Proteins and Colorectal Cancer. Current Genomics 2007,
# 8 (1), 43–49. https://doi.org/10.2174/138920207780076938.
# (13) Cao, J.; Cai, X.; Zheng, L.; Geng, L.; Shi, Z.; Pao, C. C.; Zheng, S. 
# Characterization of Colorectal-Cancer-Related cDNA Clones Obtained by
# Subtractive Hybridization Screening. Journal of Cancer Research and Clinical
# Oncology 1997, 123 (8), 447–451. https://doi.org/10.1007/BF01372549.
# (14) Bertucci, F.; Salas, S.; Eysteries, S.; Nasser, V.; Finetti, P.; 
# Ginestier, C.; Charafe-Jauffret, E.; Loriod, B.; Bachelart, L.; Montfort, J.;
# Victorero, G.; Viret, F.; Ollendorff, V.; Fert, V.; Giovaninni, M.; 
# Delpero, J.-R.; Nguyen, C.; Viens, P.; Monges, G.; Birnbaum, D.; Houlgatte, R.
# Gene Expression Profiling of Colon Cancer by DNA Microarrays and Correlation 
# with Histoclinical Parameters. Oncogene 2004, 23 (7), 1377–1391.
# https://doi.org/10.1038/sj.onc.1207262.
# (15) Bai, R.; Li, D.; Shi, Z.; Fang, X.; Ge, W.; Zheng, S. 
# Clinical Significance of Ankyrin Repeat Domain 12 Expression in Colorectal
# Cancer. Journal of Experimental & Clinical Cancer Research 2013, 32 (1), 35.
# https://doi.org/10.1186/1756-9966-32-35.

#                sessionInfo                                -----------

# R version 4.4.3 (2025-02-28)

# attached base packages:
# [1] stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
# [1] org.Hs.eg.db_3.20.0         AnnotationDbi_1.68.0        tidyr_1.3.1                
# [4] DESeq2_1.46.0               SummarizedExperiment_1.36.0 MatrixGenerics_1.18.1      
# [7] matrixStats_1.5.0           GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
# [10] IRanges_2.40.1              S4Vectors_0.44.0            sva_3.54.0                 
# [13] BiocParallel_1.40.2         genefilter_1.88.0           mgcv_1.9-3                 
# [16] nlme_3.1-168                patchwork_1.3.1             gt_1.0.0                   
# [19] vsn_3.74.0                  Biobase_2.66.0              BiocGenerics_0.52.0        
# [22] pheatmap_1.0.13             ggplot2_3.5.1               stringr_1.5.1              
# [25] dplyr_1.1.4.9000           
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.2.1        farver_2.1.2            blob_1.2.4             
# [4] Biostrings_2.74.1       fastmap_1.2.0           XML_3.99-0.18          
# [7] digest_0.6.37           lifecycle_1.0.4         survival_3.8-3         
# [10] statmod_1.5.0           KEGGREST_1.46.0         RSQLite_2.3.9          
# [13] magrittr_2.0.3          compiler_4.4.3          rlang_1.1.4            
# [16] tools_4.4.3             utf8_1.2.4              S4Arrays_1.6.0         
# [19] bit_4.5.0.1             DelayedArray_0.32.0     xml2_1.3.6             
# [22] RColorBrewer_1.1-3      abind_1.4-8             purrr_1.0.2            
# [25] withr_3.0.2             grid_4.4.3              preprocessCore_1.68.0  
# [28] fansi_1.0.6             xtable_1.8-4            colorspace_2.1-1       
# [31] edgeR_4.4.2             scales_1.3.0            cli_3.6.3              
# [34] crayon_1.5.3            generics_0.1.3          rstudioapi_0.17.1      
# [37] httr_1.4.7              DBI_1.2.3               cachem_1.1.0           
# [40] affy_1.84.0             zlibbioc_1.52.0         splines_4.4.3          
# [43] parallel_4.4.3          BiocManager_1.30.26     XVector_0.46.0         
# [46] CoprManager_0.5.7       vctrs_0.6.5             Matrix_1.7-3           
# [49] jsonlite_1.8.9          bit64_4.5.2             locfit_1.5-9.12        
# [52] limma_3.62.2            affyio_1.76.0           annotate_1.84.0        
# [55] glue_1.8.0              codetools_0.2-20        stringi_1.8.4          
# [58] gtable_0.3.6            UCSC.utils_1.2.0        munsell_0.5.1          
# [61] tibble_3.2.1            pillar_1.9.0            htmltools_0.5.8.1      
# [64] GenomeInfoDbData_1.2.13 R6_2.5.1                lattice_0.22-7         
# [67] png_0.1-8               memoise_2.0.1           Rcpp_1.0.13-1          
# [70] SparseArray_1.6.2       pkgconfig_2.0.3        








