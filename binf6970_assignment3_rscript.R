library(readxl)
library(ggplot2)
library(GGally)
library(glmnet)
library(pROC)

#loading in the data
covid_data <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")

###Data Exploration 
#looking at the structure of covid_data
str(covid_data)

#looking at the correlation matrix between the numeric features
correlation <- cor(covid_data[,c(-1,-3,-4)]) 

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
covid_data_standardized <- cbind(covid_data[,c(1,3,4)], (sapply(covid_data[,c(-1,-3,-4)], scale)))

#mean of covid_data_standardized
sapply(covid_data_standardized, mean)

#standard deviation of covid_data_standardized
sapply(covid_data_standardized, sd)


##formatting the predictors and response

x <- covid_data_standardized[,c(-1,-3)]
x$SEX <- ifelse(x$SEX == 'M', 1, 0)
x <- as.matrix(x)

y <- ifelse(covid_data_standardized$Severirty == "Mild", 0, 1)
table(y)

###Train-test
set.seed(1717)

#generating indexes for the test set
test_index <- sample(dim(x)[1],round(dim(x)[1] *.25), replace = FALSE )

#train model - 10 folds
alphas <- seq(0,1,0.01)

auc_matrix <- matrix(NA, nrow = length(alphas), ncol = 2)
colnames(auc_matrix) <- c("Training AUC", "Testing AUC")

for (i in 1:length(alphas)) {
  cvx <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = alphas[i], type.measure = "auc" )
  pred_train <- predict(cvx, newx = x[-test_index,], type = "response", s=cvx$lambda.min)[,1]
  pred_test <- predict(cvx, newx = x[test_index,], type = "response", s=cvx$lambda.min )[,1]
  auc_train <- roc(y[-test_index], pred_train)
  auc_test <- roc(y[test_index], pred_test)
  auc_matrix[i,] <- c(auc_train$auc[1], auc_test$auc[1])
}

