# Differential gene expression (DGE) analysis of RNA sequencing data

# This project is divided into three parts, each of which is in a separate file:

# • The first part is in the file “analysis_part1.R”: introduction and cleaning 
#   metadata.
# • The second part ("analysis_part2.R"): data exploration.
# • The third part ("analysis_part3.R"): differential expression analysis.


# Part 1.
#       Introduction
#       Libraries
#       Load data
#       Data Cleaning
#         Counts data
#         Metadata

#                     Introduction             --------------------------------

# Understanding and predicting cancer susceptibility is an important topic in
# medicine and biology, requiring knowledge of genes associated with disease
# development. One way to identify these genes is to compare gene expression
# data in cancer patients and healthy individuals. Some genes become more active
# (up-regulated), while others have decreased expression (down-regulated),
# leading to disruption of biological pathways [1, 2].

# The data required for such an analysis is RNA-seq data (gene expression data), 
# typically in a matrix where the rows represent genes and the columns represent 
# samples. This matrix is also linked to a metadata table containing information  
# about the experimental conditions or comparison groups. In the case of cancer 
# research, this would be gene expression data in both tumor and normal 
# (control) samples.

# The goal of this project is to explore GSE68086 dataset using unsupervised 
# machine learning algorithms and find differentially expressed genes using the
# DESEq2 package.

##                         About data  

# GSE68086 dataset is available from several Internet resources.
# ◦ Information on sample preparation, RNA processing and quality control is 
# available on the Expression Omnibus GEO platform.
#    <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68086>

# There are two GSE68086 datasets on kaggle.com;
# ◦ the first one contains counts and metadata in one file : 
# <https://www.kaggle.com/datasets/willianoliveiragibin/gene-expression-v2>

# ◦ and the second one includes 4 files with counts and metadata in txt and
# csv formats, and will be used in this analysis:
# <https://www.kaggle.com/datasets/samiraalipour/gene-expression-omnibus-geo-dataset-gse68086>

# There are also publications from ScienceDirect and Wiley on studies that used
# the GSE68086 dataset [3–5].

# A detailed description of these data can be found on the NCBI website.

# The RNA read counts obtained from the high-throughput sequencing experiment 
# included 55 samples from healthy individuals ("HC", control group) and 
# 228 samples from patients with six different malignant tumors:
#     ◦ Non-small cell lung cancer      "lung"
#     ◦ Colorectal cancer               "CRC"
#     ◦ Pancreatic cancer               "pancreas"
#     ◦ Glioblastoma                    "GBM"
#     ◦ Breast cancer                   "breast"
#     ◦ Hepatobiliary carcinomas        "hepatobiliary"   

# RNA reads are mapped to the hg19 reference genome (organism: Homo sapiens).

#                     Libraries       --------------------------------
library(dplyr)
library(stringr)

#                     Load data                --------------------------------

# The original files that were used for the project are not included in
# this repository, but available at
# https://www.kaggle.com/datasets/victorsabanzagil/polymers/data
# Files: 'GSE68086_TEP_data_matrix.csv' and 'GSE68086_series_matrix.csv'

# Both files have been saved in RDS format to reduce size and can be loaded as
# follows:
counts_matrix <- readRDS("data/counts_raw.rds")
metadata_raw <- readRDS("data/metadata_raw.rds")

#                     Data Cleaning            --------------------------------

##       Counts data                           --------------------------------

# The file 'counts_raw.rds' contains counts of RNA reads for 285 samples 
# (columns) and 57736 genes (rows of the dataset). The values in the first 
# column "X" are the ENS gene identifiers (Ensembl IDs).
dim(counts_matrix)
counts_matrix[1:5, 1:5] # raw count data

# For the analysis, we need to create a DESeq2 object, which requires a count  
# matrix with all the numeric data. So we use the 'X' column to set the row 
# names of the count matrix and remove the 'X' column.
rownames(counts_matrix) <- counts_matrix$X
counts_matrix <- counts_matrix[ , -1]

##       Metadata                              --------------------------------

# The metadata contains 46 columns with descriptions of 285 samples.  
dim(metadata_raw)

# 33 out of 46 columns in the metadata have the same value, so we will omit
# those columns. 
n_unique_values <- lapply(colnames(metadata_raw), function(i){
  # 'col_name' - column name in metadata
  # 'n_levls' - number of unique values
  data.frame(col_name = i, 
             n_levls =  length(unique(metadata_raw[ , i])) )
}) %>% bind_rows()

# number of columns that have only one unique value
sum(n_unique_values$n_levls == 1) 

# These columns contain information about samples such as: 
info_type <- n_unique_values$col_name[n_unique_values$n_levls == 1] %>% 
  str_remove("X.Sample_") %>%
  str_replace_all("_", " ")

info_type[1:10]
 
# We can look at some of these columns to get an idea of how these samples were 
# obtained and processed. 
info_values <- metadata_raw[ 1, n_unique_values$n_levls == 1] %>%
  unlist() %>% 
  str_remove_all('\\"')

# Values from columns related to data processing:
sapply(c(2, 5, 6, 10:13, 27), function(i){
  toString(c(info_type[i], info_values[i])) %>% 
    str_replace(",", ":")
})

# other columns (sample characteristics):   
colnames(metadata_raw)[n_unique_values$n_levls != 1]

# The main characteristics of the samples are given in the following columns: 
chars <- c("X.Sample_characteristics_ch1",
           "X.Sample_characteristics_ch1.1",
           "X.Sample_characteristics_ch1.2",
           "X.Sample_characteristics_ch1.3",
           "X.Sample_characteristics_ch1.4",
           "X.Sample_characteristics_ch1.5")

# Each 'chars' column has values where the type of the sample characteristics 
# is the part of the string before colon and the actual value is the rest  
# of the string
metadata_raw[1:3, chars[1:3]]

# We can extract these characteristics as part of the string before the colon. 
chars_columns <- lapply(chars, function(i){
  
  col_i = metadata_raw[ , i] # colimn i 
  
  # extract the part of the string before the colon, 
  characteristics = str_extract(col_i, "^.+:") %>%
    # remove extra characters and keep unique values
    str_remove_all('\\"|:') %>% unique()
  
  data.frame(col_name = i, # metadata column name
             # characteristic names:
             characteristics = toString(characteristics ))
}) %>% bind_rows()

# Each column has two sample characteristics: 
chars_columns

# In rows 245–246, there are mixed values, for example, 'cell type' appears  
# in the column 'X.Sample_characteristics_ch1.1' up to row 245, and in 
# 'X.Sample_characteristics_ch1.5', the values of 'cell type' appear starting 
# from row 246. The same applies to other characteristics, such as 'batch', 
# 'patient id',  'cancer type', 'tissue' and 'mutational subclass'
metadata_raw[243:248, chars]

# Let's fix this. 
# The characteristics that can be extracted from these columns are:
characteristics <- str_split(chars_columns$characteristics, ", ") %>%
  unlist() %>% 
  unique()

characteristics

# column 'X.Sample_source_name_ch1' will be used as row identifier  
length(unique(metadata_raw$X.Sample_source_name_ch1)) == nrow(metadata_raw)

# Here we create a new metadata table that will contain the sample  
# characteristics defined above in the 'characteristics' object.
suppressMessages({
  metadata_processed = lapply(1:nrow(metadata_raw), function(i){
    
    # one row in metadata with characteristics columns 
    row_i = metadata_raw[i , chars] 
    
    by_char = sapply(characteristics, function(j){
      
      # Extract value for each characteristic 
      value = str_extract(str_subset(row_i, j), ': (.+)\\"$', group = 1)  
      
    }) %>% bind_cols() %>% t()
    
    # make sure there are no spaces in the values  
    colnames(by_char) = str_replace(characteristics, " ", "_")
    rownames(by_char) = NULL
    
    # The column 'X.Sample_source_name_ch1' is used as the row identifier. 
    by_char = bind_cols(X.Sample_source_name_ch1 = metadata_raw$X.Sample_source_name_ch1[i],
                        by_char)
  }) %>% bind_rows()
})

# There is also a case inconsistency for 'cell_type' values.
metadata_processed[243:248, ]

# Therefor change all values to lowcase for variables 'tissue', 'patient_id', 
# 'cell_type' 
metadata_processed <- metadata_processed %>%
  mutate(across(c(tissue, patient_id, cell_type), tolower),
         cancer_type = ifelse(nchar(cancer_type) >=4, 
                              tolower(cancer_type), 
                              toupper(cancer_type))  )

# Before dropping columns that have only one unique value and columns with
# sample IDs, we need to set the row names for the metadata to be the same as
# the column names in the count matrix.
metadata_processed$X.Sample_source_name_ch1[8:10]
colnames(counts_matrix)[8:10]

# The names of the count data columns match the values in 
# 'X.Sample_source_name_ch1', but 'X' is added when the values start with
# numbers, and the dash is replaced with a period; there are also quotation
# marks. So we change the values in "X.Sample_source_name_ch1" to match the
# column names of the counts matrix.
metadata_processed <- metadata_processed %>%
  mutate(id = str_replace_all(str_remove_all(X.Sample_source_name_ch1, '\\"'),
                              "-", "\\.")) %>%
  mutate(id = ifelse(str_detect(id, "^\\d"), paste("X", id, sep = ""), id)) %>%
  as.data.frame()

# set row names for metadata from new column 'id'
rownames(metadata_processed) <- metadata_processed$id 

# Now let's look at the variables in the new metadata table.
# The columns 'tissue' and 'cell_type' contain only one value, describing 
# the samples that were taken from blood, thrombocytes 
unique(metadata_processed$tissue)
unique(metadata_processed$cell_type)

# The columns 'X.Sample_source_name_ch1' and 'patient_id' have unique character 
# values. 
length(unique(metadata_processed$X.Sample_source_name_ch1)) == nrow(metadata_processed)
length(unique(metadata_processed$patient_id)) == nrow(metadata_processed)

# The other variables:
unique(metadata_processed$cancer_type)
unique(metadata_processed$mutational_subclass)
unique(metadata_processed$batch)

# Also note that there are two samples with low read counts in the dataset and 
# it is recommended to exclude them from the analysis.
(descr = unique(metadata_raw$X.Sample_description))

# indices for samples with low counts:
ind_low <- which(metadata_raw$X.Sample_description == descr[2])

# these samples:
samples_low <- metadata_processed[ind_low, ] %>%
  dplyr::select( cancer_type, batch, mutational_subclass) 

samples_low

# Therefore, we remove these samples from both the metadata and the count matrix 
metadata_processed <- metadata_processed[-ind_low, ] 

counts_matrix <- counts_matrix %>%
  dplyr::select(-c(rownames(samples_low)))

# Next steps:   
#  - keep the variables 'cancer_type', 'batch' and 'mutational_subclass';
#  - change the level names to avoid spaces and other symbols like "-", "+" 
#    (and "Batch0");
#  - convert them to factors.
metadata_processed <- metadata_processed %>%
  dplyr::select(cancer_type, batch, mutational_subclass) %>%
  mutate(mutational_subclass = str_replace_all(str_remove(mutational_subclass, "\\+"), ", | ", "_" ),
         batch = str_remove(batch, "Batch0")) %>%
  mutate(across(everything(), factor))

#  - for the variable 'cancer_type' we set "HC" (healthy) as the reference level 
levels(metadata_processed$cancer_type)[1]

metadata_processed$cancer_type <- relevel(metadata_processed$cancer_type, ref = "HC")

levels(metadata_processed$cancer_type)[1]

# To create a DESeq2 object later, the sample names in the metadata and counts 
# data must be in the same order. Therefore, in the last step of data processing,
# we set the metadata row names in the same order as the column names in the 
# count matrix.
ind <- match(colnames(counts_matrix), rownames(metadata_processed))
metadata_processed <- metadata_processed[ind, ]

# And save this as metadata_processed and count_matrix:
# saveRDS(metadata_processed, "data/metadata_processed.rds")
# saveRDS(counts_matrix, "data/counts_matrix.rds")
