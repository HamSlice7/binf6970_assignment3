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

#generate fold memberships
cvxf <- cv.glmnet(x[-test_index,],y[-test_index], nfolds = 10, keep = TRUE )
table(cvxf$foldid)
foldid = cvxf$foldid

#train model - 10 folds
alphas <- seq(0,1,0.01)

auc_df <- data.frame(matrix(NA, nrow = length(alphas), ncol = 3))
colnames(auc_df) <- c("Training_AUC", "Testing_AUC", "Alpha_Values")


#making a matrix for the optimal threshold for each of the alpha-values lambdamin
sn_sp_threshold_matrix_train <- matrix(NA, nrow = length(alphas), ncol = 3)
colnames(sn_sp_threshold_matrix_train) <- c("Sensitivity", "Specificity", "Threshold")

sn_sp_threshold_matrix_test <- matrix(NA, nrow = length(alphas), ncol = 3)
colnames(sn_sp_threshold_matrix_test) <- c("Sensitivity", "Specificity", "Threshold")

#for each alpha value we get the model corresponding to lambda.min. We pick the lambda.min model as our 'best model for each alpha value. From there we can measure how well each model for each alpha value performs on the training and testing data. So there will be 101 models. For each of the 101 models we also pick the best threshold. The best performing out of the 101 models will be dubbed 'the best model' and have a corresponding threshold.


for (i in 1:length(alphas)) {
  cvx <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = alphas[i], type.measure = "auc", foldid = foldid)
  pred_train <- predict(cvx, newx = x[-test_index,], type = "response", s=cvx$lambda.min)[,1]
  pred_test <- predict(cvx, newx = x[test_index,], type = "response", s=cvx$lambda.min )[,1]
  auc_train <- roc(y[-test_index], pred_train)
  auc_test <- roc(y[test_index], pred_test)
  auc_df[i,] <- c(auc_train$auc[1], auc_test$auc[1], alphas[i])
  
  #figuring out which threshold to use via min-max approach for the training data
  sensitivity_specificity_train <- cbind(auc_train$sensitivities, auc_train$specificities)
  
  #finding the index where the min value of either sensitivity or specificity is maximized and getting the associated cutoff --> gives us the cutoff with best balance between sensitivity and specificity
  indx <- which.max(apply(sensitivity_specificity_train, 1, min))
  print(sensitivity_specificity_train[indx,])
  cutoff <- auc_train$thresholds[indx]
  sn_sp_threshold_matrix_train[i,] <- c(sensitivity_specificity_train[indx,][1],sensitivity_specificity_train[indx,][2], cutoff)
  
  #figuring out which threshold to use via min-max approach for the testing data
  sensitivity_specificity_test <- cbind(auc_test$sensitivities, auc_test$specificities)
  
  #finding the index where the min value of either sensitivity or specificity is maximized and getting the associated cutoff --> gives us the cutoff with best balance between sensitivity and specificity
  indx <- which.max(apply(sensitivity_specificity_test, 1, min))
  print(sensitivity_specificity_test[indx,])
  cutoff <- auc_train$thresholds[indx]
  sn_sp_threshold_matrix_test[i,] <- c(sensitivity_specificity_test[indx,][1],sensitivity_specificity_test[indx,][2], cutoff)
  
}


sn_sp_threshold_matrix_train[which.max(apply(sn_sp_threshold_matrix_train, 1, min)),]
sn_sp_threshold_matrix_test[which.max(apply(sn_sp_threshold_matrix_test, 1, min)),]


#finding optimal 


##confusion matrices
#training labels
table(y[-test_index])

#results from predicting on the training data
table(as.numeric(pred_train > .5))

