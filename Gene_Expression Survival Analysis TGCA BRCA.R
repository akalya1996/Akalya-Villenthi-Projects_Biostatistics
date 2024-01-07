install.packages("readxl")
library(readxl)

library(shiny)
install.packages("remotes")

remotes::install_github('rstudio/shiny')
library(remotes)
library(shiny)

#########################################
install.packages("BiocManager")
library(BiocManager)
BiocManager::install('Biostrings')
library(Biostrings)


install.packages("TCGAbiolinks")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")


library(survival)
library(survminer)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

##############################################

# getting clinical data for TCGA-BRCA cohort -------------------
# Install and load the TCGAbiolinks package

library(TCGAbiolinks)

clinical_brca <- GDCquery_clinic("TCGA-BRCA")

any(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))

which(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_brca[,c(9,39,45)]


# looking at some variables associated with survival 
table(clinical_brca$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical_brca$deceased <- ifelse(clinical_brca$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == "Alive",
                                         clinical_brca$days_to_last_follow_up,
                                         clinical_brca$days_to_death)



# get gene expression data -----------

# build a query to get gene expression data for entire cohort
query_brca_all = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")


output_brca <- getResults(query_brca_all)

# get 20 primary tissue sample barcodes
tumor <- output_brca$cases[1:500]


# # get gene expression data from 20 primary tumors 
query_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor)


# download data
GDCdownload(query_brca)


# get counts
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)

brca_matrix <- assay(tcga_brca_data, "unstranded")
brca_matrix[1:250,1:250]

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))

library(DESeq2)

# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
                              colData = coldata,
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]



# vst 
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)






# Get data for TP53 gene and add gene metadata information to it -------------
brca_tp53 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "TP53")



# get median value
median_value <- median(brca_tp53$counts)


# denote which cases have higher or lower expression than median count
brca_tp53$strata <- ifelse(brca_tp53$counts >= median_value, "HIGH", "LOW")


# Add clinical information to brca_tp53
brca_tp53$case_id <- gsub('-01.*', '', brca_tp53$case_id)
brca_tp53 <- merge(brca_tp53, clinical_brca, by.x = 'case_id', by.y = 'submitter_id')


library(knitr)

library(kableExtra)

# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)
fit
ggsurvplot(fit,
           data = brca_tp53,
           pval = T,
           risk.table = T, title = "Survival Analysis for TP53")



library(ggplot2)

# Violin Plot for TP53 gene
ggplot(brca_tp53, aes(x = strata, y = counts)) +
  geom_violin(aes(fill = strata), scale = "width", width = 0.7) +
  labs(title = "Gene Expression Distribution for TP53",
       x = "Strata",
       y = "Gene Expression") +
  theme_minimal()


#############################################################


# Get data for BRCA1 gene and add gene metadata information
brca_brca1 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "BRCA1")

# Calculate the median value for BRCA1
median_value_brca1 <- median(brca_brca1$counts)

# Denote which cases have higher or lower expression than the median count for BRCA1
brca_brca1$strata <- ifelse(brca_brca1$counts >= median_value_brca1, "HIGH", "LOW")

# Add clinical information to BRCA1 data
brca_brca1$case_id <- gsub('-01.*', '', brca_brca1$case_id)
brca_brca1 <- merge(brca_brca1, clinical_brca, by.x = 'case_id', by.y = 'submitter_id')

# Fitting survival curve for BRCA1
fit_brca1 <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_brca1)

# Plot survival curve for BRCA1 using ggsurvplot
ggsurvplot(fit_brca1,
           data = brca_brca1,
           pval = TRUE,
           risk.table = TRUE, title = "Survival Analysis for BRCA1")


#######################################################################################

# Get data for BRCA2 gene and add gene metadata information
brca_brca2 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "BRCA2")

# Calculate the median value for BRCA2
median_value_brca2 <- median(brca_brca2$counts)

# Denote which cases have higher or lower expression than the median count for BRCA2
brca_brca2$strata <- ifelse(brca_brca2$counts >= median_value_brca2, "HIGH", "LOW")

# Add clinical information to BRCA2 data
brca_brca2$case_id <- gsub('-01.*', '', brca_brca2$case_id)
brca_brca2 <- merge(brca_brca2, clinical_brca, by.x = 'case_id', by.y = 'submitter_id')

# Fitting survival curve for BRCA2
fit_brca2 <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_brca2)

# Plot survival curve for BRCA2 using ggsurvplot
ggsurvplot(fit_brca2,
           data = brca_brca2,
           pval = TRUE,
           risk.table = TRUE) +
  ggtitle("Survival Analysis for BRCA2")


####################################

library(ggplot2)
library(dplyr)

# Combine data for all three genes
all_genes_data <- bind_rows(
  mutate(brca_tp53, gene = "TP53"),
  mutate(brca_brca1, gene = "BRCA1"),
  mutate(brca_brca2, gene = "BRCA2")
)

# Violin Plot for all three genes
ggplot(all_genes_data, aes(x = strata, y = counts, fill = strata)) +
  geom_violin(scale = "width", width = 0.7) +
  labs(title = "Gene Expression Distribution for TP53, BRCA1, and BRCA2",
       x = "Strata",
       y = "Gene Expression") +
  facet_wrap(~gene, scales = "free_y") +
  theme_minimal()

