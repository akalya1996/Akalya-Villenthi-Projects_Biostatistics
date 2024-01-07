library(tidyverse)
library(survival)
library(survminer)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

# Function to perform survival analysis and generate plots
perform_survival_analysis <- function(gene_name, data_matrix, metadata, clinical_data) {
  # Get data for the specified gene and add gene metadata information
  gene_data <- data_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'gene_id') %>% 
    gather(key = 'case_id', value = 'counts', -gene_id) %>% 
    left_join(., metadata, by = "gene_id") %>% 
    filter(gene_name == gene_name)
  
  # Calculate the median value for the specified gene
  median_value <- median(gene_data$counts)
  
  # Denote which cases have higher or lower expression than the median count
  gene_data$strata <- ifelse(gene_data$counts >= median_value, "HIGH", "LOW")
  
  # Add clinical information to the gene data
  gene_data$case_id <- gsub('-01.*', '', gene_data$case_id)
  gene_data <- merge(gene_data, clinical_data, by.x = 'case_id', by.y = 'submitter_id')
  
  # Fitting survival curve
  fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = gene_data)
  
  # Plot survival curve using ggsurvplot with a title
  ggsurvplot(fit, data = gene_data, pval = TRUE, risk.table = TRUE) +
    ggtitle(paste("Survival Analysis for", gene_name))
}

# Load clinical data
clinical_brca <- GDCquery_clinic("TCGA-BRCA")

# Create an "overall survival" variable
clinical_brca$deceased <- ifelse(clinical_brca$vital_status == "Alive", FALSE, TRUE)
clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == "Alive",
                                         clinical_brca$days_to_last_follow_up,
                                         clinical_brca$days_to_death)

# Load gene expression data for the entire cohort
query_brca_all <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open"
)

output_brca <- getResults(query_brca_all)
tumor <- output_brca$cases[1:20]

query_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor
)

GDCdownload(query_brca)
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, "unstranded")

# Extract gene and sample metadata
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))

# Perform survival analysis and generate plots for TP53, BRCA1, and BRCA2
perform_survival_analysis("TP53", brca_matrix, gene_metadata, clinical_brca)
perform_survival_analysis("BRCA1", brca_matrix, gene_metadata, clinical_brca)
perform_survival_analysis("BRCA2", brca_matrix, gene_metadata, clinical_brca)
