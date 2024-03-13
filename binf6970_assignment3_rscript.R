###Data Exploration 
library(readxl)
covid_data <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")
View(covid_data)

#changing covid_data to data frame class
covid_data <- as.data.frame(covid_data)

#Looking at some summary statistics
summary(covid_data)

#viewing the dimensions of the data set
dim(covid_data)

#seeing if there are any duplicate rows
dim(unique(covid_data))

#seeing if there are any duplicate columns
dim(unique(covid_data, MARGIN = 2))

#checking to see if any of the columns have a variance of less than 1
apply(covid_data, 2, function(x) var(x) < 1)

#checking for a balanced data set
table(covid_data$Severirty)

#checking the class of each feature
apply(covid_data, 2, class)

#getting the numeric character columns and converting to numeric
covid_data_2 <- apply(covid_data[,c(2,5:32)], 2, as.numeric)

#finding the mean of the numeric columns
apply(covid_data_2, 2, mean)

#finding the standard deviations of each column
apply(covid_data_2, 2, sd)

#visually inspecting the distribution
apply(covid_data_2, 2, hist)

#z-score normalization
covid_data_scale <- apply(covid_data_2, 2, scale)

apply(covid_data_scale, 2, hist)
