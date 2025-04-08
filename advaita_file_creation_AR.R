# Breanna McBean
# April 7, 2025

# Analyze gene expression differences for conditions in the AR project and
# create file formatted for Advaita analysis.

library(tidyverse)
library(DESeq2)
library(data.table)
library(genefilter)
library(org.Hs.eg.db)
library(ggfortify)

set.seed(13)
setwd("C:/Users/bmcbe/Desktop/AR_project")

# VARIABLES TO CHANGE
chosen_line <- "MDA-MB-453"
base_condition <- "R1881"
chosen_condition <- "combo_R1881"

# FUNCTIONS

pca_plot <- function(pc_1, pc_2) {
  dds_pca <- estimateSizeFactors(dds)
  pca_matrix <- counts(dds_pca, normalized = TRUE)
  size_factors <- sizeFactors(dds_pca)
  pca_matrix <- t(t(pca_matrix) / size_factors)
  pca_matrix <- t(pca_matrix)
  sample_pca <- prcomp(pca_matrix)
  # can add label=TRUE argument to autoplot to label each point with sample
  plt <- autoplot(sample_pca, x = pc_1, y = pc_2, data = metadata,
                  colour = "Treatment", size = 3) + ggtitle(chosen_line) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 20))
  ggsave(paste0(getwd(), "/plots/", pc_1, "_", pc_2, ".png"), plot = plt)
}

# employ deseq2
run_deseq2 <- function() {
  dds <- DESeq(dds_1)
  # remove genes with on average <1 read per sample
  dds <- dds[rowSums(counts(dds)) >= dim(dds)[2]]
}

# function to create filter significantly DEG and save genes in advaita format
# NOTE: current version of advaita requires all differentially expressed genes,
# not just the significant ones for their analysis
sig_genes_and_advaita_file <- function() {
  # Order of contrast: + means numerator has more, - means denom has more
  res <- results(dds,
                 contrast = c("Treatment", chosen_condition, base_condition),
                 alpha = .05)
  res_df_tmp <- as.data.frame(res)
  res_df <- res_df_tmp[!is.na(res_df_tmp$padj), ]
  # Add entrez IDs for KEGG)
  res_df$Symbol <- mapIds(org.Hs.eg.db, key = row.names(res_df),
                          column = "SYMBOL", keytype = "ENSEMBL")
  return(res_df)
}

# END OF FUNCTIONS

# files to read
input_file <- "expression_matrix_1225-AM.tsv"
metadata_file <- "metadata_1225-AM.csv"

# read in data file and metadata file
expression_data_0 <- read_tsv(input_file)
metadata <- read.csv(metadata_file, sep = ",")

# cell lines, regimens, and batches in experiment
cell_lines <- unique(metadata$Cell_line)
treatments <- unique(metadata$Treatment)

# switch gene_IDs column to row names
expression_data <- expression_data_0 %>% column_to_rownames(var = "gene_ID")
all_genes <- row.names(expression_data)

# create factors for cell line and regimen
metadata$Cell_line <- factor(metadata$Cell_line, cell_lines)
metadata$Treatment <- factor(metadata$Treatment, treatments)
samples_meta <- distinct(metadata, metadata$Sample_ID, .keep_all = TRUE)
metadata <- samples_meta

# split dataframe for chosen cell line/regimen and by treatment
metadata_tmp <- samples_meta[samples_meta$Cell_line == chosen_line, ]
metadata <- metadata_tmp[metadata_tmp$Treatment %in% c(chosen_condition,
                                                       base_condition), ]
expression_data <- expression_data[, unique(metadata$Sample_ID)]
treatments <- c(base_condition, chosen_condition)
metadata$Treatment <- factor(metadata$Treatment, treatments)

# create deseq2 object
dds_1 <- DESeqDataSetFromMatrix(countData = expression_data,
                                colData = metadata,
                                design = ~Treatment)

# call function to run deseq2
dds <- run_deseq2()

pca_plot(1, 2)

# call function to create filter significantly DEG and save genes in advaita-
# friendly format
sig_res_df <- sig_genes_and_advaita_file()
