###Data Exploration 
library(readxl)
covid_data <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")

#looking at the structure of covid_data
str(covid_data)


#checking the class of each feature
sapply(covid_data, class)

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

#finding the mean of the numeric columns
sapply(covid_data, mean)

#finding the standard deviations of each column
sapply(covid_data, sd)

#visually inspecting the distribution
sapply(covid_data, function(x) {if (class(x) == "numeric"){hist(x)}})


#calculating the means for mild and severe groups exclusively 

covid_data_mild <- covid_data[covid_data$Severirty == "Mild",]
covid_data_severe <- covid_data[covid_data$Severirty == "Severe",]

#calculating mean for covid_data_mild
mean_mild <- as.data.frame(sapply(covid_data_mild, mean))

#calculating mean for covid_data_severe
mean_severe <- as.data.frame(sapply(covid_data_severe, mean))

#looking at the differences in mean
abs(mean_mild - mean_severe)

#z-score normalization
covid_data_scale <- apply(covid_data_2, 2, scale)

apply(covid_data_scale, 2, hist)

