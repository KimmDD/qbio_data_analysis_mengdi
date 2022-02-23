# Exercise 1.1
BiocManager::install("DESeq2")
# Load packages here
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
# set our working directory
setwd("~/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi")
# create our sum_exp object again
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method


sum_exp <- GDCprepare(query)


# Exercise 1.2
# code goes here!

# 1. Identify which patients have NA in their age as before
bool_age_na = is.na(sum_exp$age_at_diagnosis)
bool_age_na
# 2. Make a copy of the clinical and counts data; fill in here
patients_data = colData(sum_exp)  # contains the clinical data
counts = assays(sum_exp)$"HTSeq - Counts" # contains the counts data
# 3. Convert counts into a dataframe; fill in here
# 4. Remove the NA patients from both patients_data and counts
patients_data = patients_data[!bool_age_na, ]
patients_data$age_at_diagnosis = patients_data$age_at_diagnosis / 360
patients_data$age_at_diagnosis
counts = counts[ , !bool_age_na]
# 5. Create the age_category column
patients_data$age_category = ifelse(patients_data$age_at_diagnosis < 50, "young", "old")
patients_data$age_category
# 6. Turn the age_category column into a factor
patients_data$age_category = factor(patients_data$age_category, levels = c("young", "old"))
patients_data$age_category
# We need to specify the factor order to make sure we compare young over old
# If you capitalized young and old differently, maske the levels=c(...) parameter vector match



# Exercise 1.3
#1. Convert the rownames of counts to make them human-readable
# make sure that your Ensembl row names match perfectly to the data in rowRanges before replacing
if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}

counts
# 2. Compute the sum across the row of the counts data frame
counts_row_sums = rowSums(counts)
counts_row_sums
# 3. Identify the genes that have fewer than 10 reads
low_counts_mask = counts_row_sums < 10
sum(low_counts_mask)
# 4. Remove these lowly expressed genes from counts
counts = counts[!low_counts_mask, ]

is.na(patients_data)
ncol(counts)
nrow(counts)
nrow(patients_data)
patients_data$age_category

# Exercise 2.1
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patients_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

results

# Exercise 2.2
# look at the results object
str(results)
head(results)


# Exercise 2.3
# Create a new vector, row_order
row_order = order(results$padj)
# Looking at the previous example, use row_order to sort results
results = results[row_order,]
# Then, look at the first 20 rows of results using head()
head(results,20)
# Pick a significantly differentially expressed gene
# ENSG00000252010, it is more highly expressed in old patients
# Small Cajal Body-Specific RNA 5, the function is bilirubin measurement, calcium measurement, etc.

results$log2FoldChange

# Exercise 2.4
log2FoldChange_threshold = 1
padj_threshold = 0.05
# Identify the genes that have a log2FoldChange that is greater than your log2FoldChange_threshold
genes_greater_log2 = results[results$log2FoldChange > log2FoldChange_threshold, ]
genes_greater_log2
# Identify the genes that have a padj that is less than your pad_threshold
genes_smaller_padj = results[results$padj > padj_threshold, ]
genes_smaller_padj 
# Subset the results data frame, selecting the genes that satisfy BOTH criteria
subset.data.frame(results, genes_greater_log2, genes_smaller_padj)
