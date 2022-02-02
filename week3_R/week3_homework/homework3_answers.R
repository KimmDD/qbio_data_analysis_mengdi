# Exercise 1.1 

# access the attenu data set 
attenu

# Identify the rows in which there was no reported station in the station column
# The answer is 79  81  94  96  99 107 108 114 116 118 123 126 128 155 156 160
which(is.na(attenu$station))

# Create a copy of attenu called attenu_cleaned
# where rows that are missing station information are not included.
attenu_cleaned <- na.omit(attenu)

#Print the first 6 rows of attenu_cleaned (using head()) 
# and its dimensions (using dim())
head(attenu_cleaned,6)
dim(attenu_cleaned)



# Exercise 1.2

# Make a copy of the Theoph data set called Theoph_2
Theoph_2 <- Theoph

# Identify what the median treatment dose is
str(Theoph_2)
median_dose = median(Theoph_2$Dose)

# using ifelse(), add a new column called Dose_Class to Theoph_2
#where the dose is classified as high 
# if it is above or equal to the median dose, and 'low otherwise.
Theoph_2$Dose_Class <- NA
Theoph_2$Dose_Class = ifelse(Theoph_2$Dose < median_dose, Theoph_2$Dose_Class <- "low",
       Theoph_2$Dose_Class <- "high")

# print the first 6 rows and its dimensions.
head(Theoph_2,6)
dim(Theoph_2)



# Exercise 1.3

# Read in the Starbucks nutrition data 
# and store it into variable called "starbucks"
starbucks <- read.csv("/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week3_R/week3_homework/starbucks.csv", header = TRUE)

# Determine which values contain NA using is.na()
!is.na(starbucks)

# Create a boolean vector, is_row_empty 
# Sorry I do not know how to the rest of parts



# Exercise 1.4

# Read the data in using the read.csv() function.
Batting <- read.csv("/Users/chaimengdi/Desktop/QBIO/QBIO490/qbio_data_analysis_mengdi/week3_R/week3_homework/Batting.csv", header = TRUE)

#  Identify the number of players who scored three or more home runs in a given year.
Batting_subset = Batting[Batting$HR >= 3, ]
nrow(Batting_subset)

# Plot the number of homeruns vs. year
plot(Batting$yearID, Batting$HR)

# Create a new data frame, containing players from the LA Angels.
Batting_LA <- Batting[Batting$teamID == 'LAA', ]
plot(Batting_LA$yearID, Batting_LA$HR)

# Repeat step 4, except subset the original batting data to include players from either “ATL” or “PIT”
Batting_sub <- Batting[Batting$teamID == 'ATL' | Batting$teamID == 'PIT' , ]
plot(Batting_sub$yearID, Batting_sub$HR)




# Exercise 2.1
#  Iris contains three plant species and four features measured for each sample. 
# These quantify the morphologic variation of the iris flower 
# in its three species, all measurements given in centimeters.
# it contains 150 observations
# The feature is about Sepal.length, Sepal.width, Petal.length, Petal.width
# and Species 

# Exercise 2.2
# I think Sepal.length, Sepal.width, Petal.length, Petal.width are cintunuous
# variables, and Species are categorical. Sepal.length, Sepal.width, 
# Petal.length, Petal.width are all integers, and Species are strings.

# Exercise 2.3
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)

# interesting: for Sepal's length and width, it follows the normal distribution
# However, the Petal'length and width has more on the two side and less at the center

# Exercise 2.4

# Find the mean sepal width in the data frame, and store it in a variable
mean_sepal_width = mean(iris$Sepal.Width)
# make a copy of iris called iris_copy
iris_copy <- iris
# Create a new vector using the ifelse() function to fill out values, comparing Sepal.Width to the mean value
# Create a new column in the iris data set (using the dollar sign) and assign the new vector to that column
iris_copy$Width_Species = ifelse(iris_copy$Sepal.Width < mean_sepal_width, "narrow-sepaled", "wide-sepaled")


# Create a boxplot plotting the Sepal.Width based on the new column
boxplot(iris_copy$Width_Species ~ iris_copy$Sepal.Width)







 