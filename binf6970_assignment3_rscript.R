###Data Exploration 
library(readxl)
covid_data <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")
View(covid_data)

#changing covid_data to data frame class
covid_data <- as.data.frame(covid_data)

#Looking at some summary statistics
summary(covid_data)

