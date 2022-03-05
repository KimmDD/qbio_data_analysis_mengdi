# General Concepts

# What is TCGA and why is it important?
# TCGA is The Cancer Genome Atlas. 
# It is organized by the National Cancer Institute (NCI) and National Human Genome Research Institute (NHGRI)
# It cataloged data from over 20,000 samples of 33 cancer types
# These data is the best source for scientists to research on the relationship between genes and cancer 

# What are some strengths and weaknesses of TCGA?
# strength: contains many data such as clinical data, genomic mutation data, mRNA transcriptomic data and so on 
# weakness: the sample size is a little bit small, leading to some analysis may be not accurate enough

# How does the central dogma of biology (DNA → RNA → protein) relate to the data we are exploring?
# the dogma of biology relates to our data in the database, for example the mRNA transcriptomic data. 
# The mRNA transcriptomic data we use could draw the plot connecting the RNA and protein together
# And find something that we miss.


# Coding Skills
# What commands are used to save a file to your GitHub repository?
# git commit 

# What command must be run in order to use a package in R?
# library


# What is boolean indexing? What are some applications of it?
# the boolean indexing returns true or false and we could use this to indicate the NA value and exclude the NA values in database
# Also we could select certain data based on given row or column

# Draw out a dataframe of your choice. Show an example of the following and explain what each line of code does. 
library(TCGAbiolinks)
# show the column name of dataframe clinic
colnames(clinic)
# pick age_at_initial_pathologic_diagnosis column 
clinic$age_at_initial_pathologic_diagnosis
# the ifelse statement creates a new column that if the patient's age is less than 50,
# they will be identified as young, otherwise, they will be identified as old
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "young", "old")
# This boolean indexing examines the patients' race is empty or not. If the race of patients is NA,
# it will show true. Otherwise, it will show false 
is.na(clinic$race_list)



