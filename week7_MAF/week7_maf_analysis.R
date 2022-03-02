# Exercise 1.1
# install maftools, load packages, then
BiocManager::install("maftools")
library(TCGAbiolinks)
library(maftools)
# set our working directory
setwd("~/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi")


# Exercise 1.2
# replace the path!
# fread is a way to read in a data frame, but it's much faster than the built-in read.csv()
# by default, it reads it in as a slightly different data type than a data frame, so we set the data.table flag to false
clinic <- read.csv("~/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/coad_clinical_data.csv", row.names=NULL)

# rename the patient barcode to make it work with maftools
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"



# Exercise 1.3
colnames(clinic) # 76 long
colnames(clinic) == "bcr_patient_barcode" # 76 long, contains boolean, no true value



# Exercise 1.4
# fill in!!
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)



# Exercise 1.5
# 2. Locate the appropriate csv file
setwd("/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi")
list.files("GDCdata")

# 3. Use fread to read in the *mutation* data you just downloaded (refer to the fread syntax from the previous example)
maf_dataframe <- read.csv("~/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", row.names=NULL)

# 4. Read in your clinic data frame
# you also need to rename the bcr_patient_barcode again
# refer to the previous example!
clinic <- read.csv("~/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/coad_clinical_data.csv", row.names=NULL)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

# re-create the maf_object as before
maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)



# Exercise 2.1
# Just type the name of the MAF object and run that line alone
maf_object
# Call str() on it
str(maf_object)
# Access the data and clinical.data data frames stored in MAF object
maf_object@data
colnames(maf_object@data)
maf_object@clinical.data
colnames(maf_object@clinical.data)



# Exercise 3.1
# play with the number of genes with the "top" argument
oncoplot(maf = maf_object,
         top = 10) 

oncoplot(maf = maf_object,
         top = 5) 


# Exercise 3.2
# The gene we choose is APC 
# APC protein acts as a tumor suppressor
# which means that it keeps cells from growing and dividing too fast or in an uncontrolled way
# It helps control how often a cell divides


# Exercise 3.3
# 1. Write to clinic again
clinic = maf_object@clinical.data

# 2. Create the young_patients_ids vector
clinic$Tumor_Sample_Barcode
young_patients_ids = clinic$Tumor_Sample_Barcode[clinic$age_at_initial_pathologic_diagnosis < 50]

# 3. Create another 
young_maf = subsetMaf(maf = maf_object,
                      tsb = young_patients_ids)

# 4. Repeat steps 2-3 to create an old_maf! Can you do it in one line?
old_patients_ids = clinic$Tumor_Sample_Barcode[clinic$age_at_initial_pathologic_diagnosis >= 50]
old_maf = subsetMaf(maf = maf_object,
                      tsb = old_patients_ids)


# Exercise 3.4
coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "young", 
           m2Name = "old")




# Exercise 3.5
# The gene mutation number is mixed in old and young patients and that makes sense
# because the young and old patients have different genomic expression related to aging 
# so the different types of mutation may more occur in different groups of people 


# Exercise 3.6
# pick a gene to look at!
lollipopPlot(maf_object, gene = "APC")


# Exercise 3.7
# pick a gene to look at!
lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name= "young_patients",
              m2_name = "old_patients",
              gene = "APC")
# The gene is more commonly mutated in old patients
# There are more mutation in a specific region of the protein
# That is because the protein consist with amino acid sequence 
# and may some amino acid sequence causes this gene mutation easily
# The Missense_mutation, Splice_site, and Nonsense_mutation are more common
# The distribution of missense_mutation is at the end of the protein, and
# the distribution of nonsense_mutation and splice_site are more often at the beginning of the protein


# Exercise 3.8
# There are six samples that have both a mutation in gene A and gene B


# Exercise 3.9
b = 7
c= 2
d = 35
e= 37
f = 42

# Exercise 3.10
# geneA is TP53
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

# geneB is KRAS
geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")


# Exercise 3.11
geneA_maf
geneB_maf # subsetMaf() examines the class of MAF
# each sample may have more than one mutation 
maf_object@clinical.data # The number of sample is different
# I guess it is because not every sample contains the gene of TP53 or KRAS?



# Exercise 3.12
# 1. Access the barcodes of the patients with mutations in genes A and B
# bc stands for barcode
mut_bc_geneA = geneA_maf@clinical.data$Tumor_Sample_Barcode
mut_bc_geneB = geneB_maf@clinical.data$Tumor_Sample_Barcode

# 2. Get the lengths of these two vectors
num_mut_geneA = length(mut_bc_geneA)
num_mut_geneB = length(mut_bc_geneB)

# 3. Fill in the intersect here! Then get the nubmer of patients
mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)



# Exercise 3.13
num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
num_mut_geneB_only = num_mut_geneB - num_mut_geneAB



# Exercise 3.14
maf_object #total number is 399
num_mut_neither = 399 - num_mut_geneB - num_mut_geneA_only

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_mut_neither), 
                       nrow=2)

# view the contingency table
contig_table



# Exercise 3.15
fe_results <- fisher.test(contig_table)
fe_results

