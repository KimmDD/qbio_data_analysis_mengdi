if(!require(BiocManager))install.packages("BiocManager")

if(!require(TCGAbiolinks)) {
  BiocManager::install("TCGAbiolinks")
}

library(TCGAbiolinks)

# You will use variations of this command for every data type
clin_query <- GDCquery(project = "TCGA-COAD",
                       data.category = "Clinical",
                       file.type = "xml")


# ONLY RUN THIS ONCE
# GDCdownload(clin_query)

clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")

# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"


# Exercise 1.1
# look at the clinic object here
str(clinic)
head(clinic)

# Exercise 1.2
# Look at the clinic object using colnames
colnames(clinic)
# # Access any column using the $ syntax, for example clinic$vital_status
clinic$vital_status

# Exercise 2.1
plot(clinic$age_at_initial_pathologic_diagnosis,clinic$weight, xlab="age", ylab="weight")

# Exercise 2.2
unique(clinic$race_list) # This shows you all the UNIQUE entries of the race column
#the mar argument in par() sets the plot margins. 
# As the race names may be long, we want to have a large bottom margin. 
# The syntax is par(mar = c(bottom, left, top, right)), where mar is margins
# How does the plot change if you change these margins? (shift to left or right)
par(mar=c(10,1,1,1))

# there's a few extra parameters to make the plot display better
boxplot(clinic$age_at_initial_pathologic_diagnosis ~ clinic$race_list,  # fill in variables to plot here
        las = 2, 
        cex.axis = 0.5)

# Exercise 2.3
# Determine the number of patients without race information
sum(is.na(clinic$race_list))
# replace these "" with "No data"
clinic$race_list = as.character(clinic$race_list)
clinic$race_list[is.na(clinic$race_list)] <- "No data"

# Exercise 2.4
# find the min, max, mean, and median ages
min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)

# also, run the summary() function
summary(clinic$age_at_initial_pathologic_diagnosis)

# Exercise 2.5
sum(clinic$age_at_initial_pathologic_diagnosis < 50)
sum(clinic$age_at_initial_pathologic_diagnosis >= 50)

# Exercise 2.6
young_patient_ids = which(clinic$age_at_initial_pathologic_diagnosis < 50)
old_patient_ids = which(clinic$age_at_initial_pathologic_diagnosis >= 50)

# Exercise 2.7
# create the new column
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis <50, "young", "old")

# Exercise 2.8
clinic[1,1] #This is the top left entry of the dataframe. R has "one-based indexing"
clinic[1,] # This is the first row of the dataframe; the blankspace after the comma is to contain all columns in row one
clinic[2:5,] # This is the second to fifth rwo of the dataframe
clinic[,3] # the is the third column of the dataframe; contains all rows from column 3

# Exercise 2.9
# create two new data frames, young_clinic and old_clinic
# hint: use the syntax from before
split_clinic <- split(clinic, f = clinic$age_category)

old_clinic = split_clinic$old
young_clinic = split_clinic$young

# Exercise 2.10
young_clinic_one_line <- clinic[clinic$age_at_initial_pathologic_diagnosis < 50,]
young_clinic_one_line

# Check if your dimensions are the same! dim() gives dimensions of a data frame in number_of_rows, number_of_cols
identical(dim(young_clinic), dim(young_clinic_one_line))



# install and load your packages in here!
install.packages("survival")
install.packages("survminer")

library(survival)
library(survminer)

# Exercise 3.1
# clean the NA values here!
clinic$days_to_death = ifelse(!is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death)

# Exercise 3.2
# create the death_event column here
clinic$death_event = ifelse(clinic$vital_status == "Alive", 0, 1)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinic$race_list, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

# save the plot as a png
# if you want to present the data, change the strata labels!
ggsave("/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/kmplot_by_race.png", plot = p, width = 12, height = 9)

# Exercise 3.3
# it is decreasing and similar to stairs
# Question: why at some points, the survival probability decreases dramatiacally, whereas at some points is not?
# I think it should add the units for x-axis and y-axis.

# Exercise 4.1
# change the file path! make sure it's in your week4 folder
# we set row.names to false since the rows don't contain any relevant info
write.csv(clinic, "/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/coad_clinical_data.csv", row.names = F)



