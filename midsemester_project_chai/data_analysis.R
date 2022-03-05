library(TCGAbiolinks)
library(SummarizedExperiment)
setwd("~/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method

sum_exp <- GDCprepare(query)

# Question 1: 
# What populations (men/women, young/old, etc.) tend to over gene CD99? 

# Find the Ensembl gene ID of gene CD99
gene_CD99_mask =  (rowData(sum_exp)$external_gene_name == "CD99")
ensembl_Cd99 = rowData(sum_exp)$"ensembl_gene_id"[gene_CD99_mask]

# Find the min and max number of counts of CD99
summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_Cd99, ]  )

# find the NA values on gender
bool_gender_na = is.na(colData(sum_exp)$gender)
# exclude the NA values in gender
gender_cat_no_NAs = sum_exp$gender[!bool_gender_na]
# Count the genes with gender
gene_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl_Cd99,!bool_gender_na]

# Create the boxplot with gender on x-axis and the counts of gene CD99 on y-axis 
boxplot(gene_counts ~ gender_cat_no_NAs, 
        main = "CD99 counts on male vs. female",
        xlab = "Gender", 
        ylab = "Gene Count",
        cex.axis = 1,
        col="orange")


# Question 2:
# How does gene CD99 affect the survival time? 
library(survival)
library(survminer)
sum_exp_dataframe<- as.data.frame(colData(sum_exp))


# replace the NA balue in days_to_death to days_to_last_follow_up
sum_exp_dataframe$days_to_death = ifelse(is.na(sum_exp_dataframe$days_to_death), 
                                       sum_exp_dataframe$days_to_last_follow_up, sum_exp_dataframe$days_to_death)
# make the days_to_death column numeric 
sum_exp_dataframe$days_to_death = as.numeric(sum_exp_dataframe$days_to_death)
# create the death_event column here
sum_exp_dataframe$death_event = ifelse(sum_exp_dataframe$vital_status == "Alive", 0, 1)
# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = sum_exp_dataframe$days_to_death, 
                    event = sum_exp_dataframe$death_event)

# We then create a fit object to do the gender analysis 
gender_fit <- surv_fit( surv_object ~ sum_exp_dataframe$gender, data = sum_exp_dataframe )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# draw the final plot
p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

# create the new category with low expression of CD99 patients and high expression patients
sum_exp_dataframe$expression_category = ifelse(assays(sum_exp)$"HTSeq - Counts"[ensembl_Cd99, ]< 6956, "low", "high")

# create the expression fit
expression_fit <- surv_fit( surv_object ~ sum_exp_dataframe$expression_category, data = sum_exp_dataframe )
# draw the survive plot comparing low expression and high expression 
survplot2 = ggsurvplot(expression_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# draw the final plot 
m = survplot2$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

m
# create the survive plot for male with high expression and low expression 
male_group = sum_exp_dataframe[sum_exp_dataframe$gender == "male", ]
# make the days_to_death column numeric for male
male_group$days_to_death = as.numeric(male_group$days_to_death)
# create the death_event column here for male 
male_group$death_event = ifelse(male_group$vital_status == "Alive", 0, 1)
# We initialize a 'survival' object first, which contains the data we need.
male_surv_object <- Surv(time = male_group$days_to_death, 
                    event = male_group$death_event)

male_fit <- surv_fit( male_surv_object ~ male_group$expression_category, data = male_group )

male_survplot = ggsurvplot(male_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# draw the final plot
male = male_survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
male

# create the survive plot for female with high expression and low expression 
female_group = sum_exp_dataframe[sum_exp_dataframe$gender == "female", ]
# make the days_to_death column numeric for male
female_group$days_to_death = as.numeric(female_group$days_to_death)
# create the death_event column here for male 
female_group$death_event = ifelse(female_group$vital_status == "Alive", 0, 1)
# We initialize a 'survival' object first, which contains the data we need.
female_surv_object <- Surv(time = female_group$days_to_death, 
                         event = female_group$death_event)

female_group = sum_exp_dataframe[sum_exp_dataframe$gender == "female", ]

female_fit <- surv_fit(female_surv_object ~ female_group$expression_category, data = female_group )

female_survplot = ggsurvplot(female_fit, 
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                           legend = "right")
# draw the final plot
female = female_survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
female


# Question 3:
# How's the CD99 mutation related to sex group 
library(maftools)
maf_dataframe <- read.csv("~/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", row.names=NULL)
clinic <- read.csv("~/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/coad_clinical_data.csv", row.names=NULL)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
maf_object <- read.maf(maf = maf_dataframe, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

# Write to clinic
clinic = maf_object@clinical.data
# create the male vector
male_patients_ids = clinic$Tumor_Sample_Barcode[clinic$gender == "MALE"]
# create the female vector
female_patients_ids = clinic$Tumor_Sample_Barcode[clinic$gender == "FEMALE"]
# create the male submaf
male_maf = subsetMaf(maf = maf_object,
                      tsb = male_patients_ids)
# create the female submaf
female_maf = subsetMaf(maf = maf_object,
                     tsb = female_patients_ids)
# draw the lollipopplot for CD99 gene
lollipopPlot(maf_object, gene = "CD99")
# draw the lollipopplot for CD99 gene comparing with male and female
lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name= "male_patients",
              m2_name = "female_patients",
              gene = "CD99")

# We could draw the conclusion that the gene CD99 mutation occurs in male instead of female 






