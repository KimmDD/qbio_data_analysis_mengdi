# Exercise 1.1
# # Install the package here
BiocManager::install("SummarizedExperiment")
# Load packages here
library(TCGAbiolinks)
library(SummarizedExperiment)
# set your working directory here
setwd("~/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi")


# Exercise 2.1
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method

# only need to download the data once! Comment this out once you have completed it once
GDCdownload(query)

sum_exp <- GDCprepare(query)

# Exercise 2.2
assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]

# Exercise 2.3
# the number of row for rowData(sum_exp) and assays(sum_exp)$"HTSeq - Counts" is same
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")


# Exercise 2.4
str(colData(sum_exp))
head(colData(sum_exp))

# Exercise 2.5
str(colnames(sum_exp))
head(colnames(sum_exp))
# the name of the column with the age of the patients
sum_exp$age_at_index

# Exercise 2.6
# unit of column is day(s)
sum_exp$age_at_diagnosis[1:10]

# Exercise 2.7
colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis 
#finish this line using a math operator to convert to years
sum_exp$age_at_diagnosis = sum_exp$age_at_diagnosis / 360

# Exercise 2.8
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis < 50, "Young", "Old")


# Exercise 2.9
head(rowData(sum_exp))
dim(rowData(sum_exp))


# Exercise 2.10
# Try with your gene name and the external_gene_name column!
# Hint: look at the previous example
"MSH2" %in% rowData(sum_exp)$external_gene_name
"MSH6" %in% rowData(sum_exp)$external_gene_name



# Exercise 2.11
# Look at rows here!
assays(sum_exp)$"HTSeq - Counts"[30:35, 20:25]


# Exercise 2.12
#STEP 1
geneA_id_mask = (rowData(sum_exp)$external_gene_name == "MSH2") #what should go before the $ in this line?
sum(geneA_id_mask) #BEFORE running this line, guess the result. 
# Hint: How many TRUE's should be in this mask? How many times should your gene name appear?

#STEP 2
ensembl_geneA = rowData(sum_exp)$"ensembl_gene_id"[which(geneA_id_mask)] #fill in the dotted lines. 

#STEP 3: Repeat for geneB
geneB_id_mask = (rowData(sum_exp)$external_gene_name == "MSH6")
ensembl_geneB = rowData(sum_exp)$"ensembl_gene_id"[which(geneB_id_mask)]


# Exercise 2.13
# It is the column

# Exercise 2.14
min(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ])  #On which side of the comma does ensembl_geneA go?
# find max of geneA
max(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ])
summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ]  ) # Use the same thing as for the min() function, but for your second gene.


# Exercise 2.15
plot(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ],
     assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ],
     xlab = "Gene A", # remember to rename axes!
     ylab = "Gene B"
)


# Exercise 2.16
bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na


# Exercise 2.17
age_cat_no_NAs = sum_exp$age_category[!bool_age_na]
age_cat_no_NAs


# Exercise 2.18
# #use the length() function here
length(age_cat_no_NAs)
#create your math equation here. Hint: You can either add two of the values to equal the third, or subtract two values to equal the third. 
dim( colData(sum_exp) ) #gives number of rows then number of columns
num_na + length(age_cat_no_NAs) == 521 #the number of rows in colData(sum_exp)


# Exercise 2.19
# determine the number of patients here
dim(assays(sum_exp)$"HTSeq - Counts") # there are 521 patients 

# Exercise 2.20
identical( rownames(colData(sum_exp)),colnames(assays(sum_exp)$"HTSeq - Counts"))
# Fill in the dotted lines with "row" or "col"

gene_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,!bool_age_na]

length(gene_counts)
# Exercise 2.21
# Check here using the `length()` function and `==`
length(age_cat_no_NAs) == length(gene_counts)


# Exercise 2.22
# make your boxplot here
# please separate your arguments onto separate lines to make it easier to read!
boxplot(gene_counts ~ age_cat_no_NAs, 
        xlab = "Age Category", 
        ylab = "Gene Count")


# Exercise 3.1
# 1
assays(sum_exp)$"HTSeq - Counts" # the row represents the genes, column represents the alternative gene names
# 2
# first we could get rid of the NA values, then tracking the index of the row 
# then find it in the counts data frame
# the genes information are stored 
# This data frame is included into the counts data frame
# 3
# first we could get rid of the NA values, then tracking the index of the row 
# then find it in the counts data frame
# the anternative names of genes are stored 
# This data frame is included into the counts data frame





