clinic <- read.csv("/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/coad_clinical_data.csv")

# Written Activity

# 1
# Categorical variables are variables that classify observations into groups.
# Example: countries, years, genders. etc.

# discrete variables are variables for which the values 
# it can take are countable and have a finite number of possibilities
# Example: the number of feedbacks you received 

# Continuous variables are numeric variables that 
# have an infinite number of values between any two values
# Example: the length of a part of something 

# 2
clinic_age = clinic$age_at_initial_pathologic_diagnosis

# 3
# The age is from the participants and collected at the first diagnosis
# It is a continuous variable 

# 4
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544764/
# Summary: This article challenges the idea that cancer cannot be prevented 
# among older adults by examining different aspects of the relationship 
# between age and cancer. 
# https://www.nature.com/articles/s41416-019-0721-1
# Summary: This article describe the mechanisms that we and other animals 
# have evolved to avoid cancer through reproductive years, 
# including by limiting the accumulation of mutations, particularly if 
# these mutations alter cellular phenotypes, and by maintaining youthful tissue 
# landscapes so as to favour the status quo. 

# 5
clinic_height = clinic$height
# The weight describes the participants' weights
# and is collected when the participants first have diagnosis
# It is a continuous variable 

# 6
# the age and gender are randomly assigned and has no closed pattern or distribution 
# People who have older age may have lower survival rate in colorectal cancer
# Male has higher risk in survival in coloreactal cancer when female 

# 7
# When I finished analyze both age and weight, I found that the graphs are very 
# different from the kmplot_by_race. I think the reason is that for age and weight,
# they are only have two categories. However, for race-list, it contains many kinds of races.


# Coding Activity
# 1
clinic$age_categorical = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "young", "old")
boxplot(clinic$height ~ clinic$age_categorical, xlab="age", ylab="weight")
clinic$weight_categorical = ifelse(clinic$weight < 150, "low weight", "high weight")

# 2
library(survival)
library(survminer)

# clean the NA values here!
clinic$days_to_death = ifelse(!is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death)
# create the death_event column here
clinic$death_event = ifelse(clinic$vital_status == "Alive", 0, 1)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
age_fit <- surv_fit( surv_object ~ clinic$age_categorical, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(age_fit, 
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
ggsave("/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/kmplot_by_age.png", plot = p, width = 12, height = 9)


# 3
# We then create a fit object
weight_fit <- surv_fit( surv_object ~ clinic$weight_categorical, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(weight_fit, 
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
ggsave("/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week4_clinical/kmplot_by_weight.png", plot = p, width = 12, height = 9)


