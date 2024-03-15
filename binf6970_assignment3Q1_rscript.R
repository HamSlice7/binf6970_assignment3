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

#looking at sex balance
table(covid_data$SEX)

#looking at data balance
table(covid_data$Severirty)

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
x$SEX <- as.factor(x$SEX)
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
sn_sp_threshold_matrix_train <- matrix(NA, nrow = length(alphas), ncol = 4)
colnames(sn_sp_threshold_matrix_train) <- c("Sensitivity", "Specificity", "Threshold", "Alpha_Value")

sn_sp_threshold_matrix_test <- matrix(NA, nrow = length(alphas), ncol = 4)
colnames(sn_sp_threshold_matrix_test) <- c("Sensitivity", "Specificity", "Threshold", "Alpha_Value")

#for each alpha value we get the model corresponding to lambda.min. We pick the lambda.min model as our 'best model for each alpha value. From there we can measure how well each model for each alpha value performs on the training and testing data. So there will be 101 models. For each of the 101 models we also pick the best threshold. The best performing out of the 101 models will be dubbed 'the best model' and have a corresponding threshold.


for (i in 1:length(alphas)) {
  #training the model
  cvx <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = alphas[i], type.measure = "deviance", foldid = foldid)
  #using the trained model to make predictions of the training data
  pred_train <- predict(cvx, newx = x[-test_index,], type = "response", s=cvx$lambda.min)[,1]
  #using the trained model to make predictions for the testing data 
  pred_test <- predict(cvx, newx = x[test_index,], type = "response", s=cvx$lambda.min )[,1]
  #calculating the area under the curve for the training predictions 
  auc_train <- roc(y[-test_index], pred_train)
  #calculating the area under the curve for the testing predictions 
  auc_test <- roc(y[test_index], pred_test)
  #appending the area under of the curve for each model to auc_df
  auc_df[i,] <- c(auc_train$auc[1], auc_test$auc[1], alphas[i])
  
  ###training data
  
  #finding the index associated with the threshold that yields the most balanced sensitivity-specificity via the min-max approach (where the min value of either sensitivity or specificity is maximized) for the training data. Each threshold has an associated sensitivity and specificity.
  sensitivity_specificity_train <- cbind(auc_train$sensitivities, auc_train$specificities)
  indx <- which.max(apply(sensitivity_specificity_train, 1, min))
  #assigning the identified threshold from min-max to the cutoff variable
  cutoff <- auc_train$thresholds[indx]
  #assign the sensitivity, specificity, and the cutoff for each for each identified model of every alpha value to sn_sp_threshold_matrix_train
  sn_sp_threshold_matrix_train[i,] <- c(sensitivity_specificity_train[indx,][1],sensitivity_specificity_train[indx,][2], cutoff, alphas[i])
  
  ###testing data
  
  ##finding the index associated with the threshold that yields the most balanced sensitivity-specificity via the min-max approach (where the min value of either sensitivity or specificity is maximized) for the testing data. Each threshold has an associated sensitivity and specificity.
  sensitivity_specificity_test <- cbind(auc_test$sensitivities, auc_test$specificities)
  indx <- which.max(apply(sensitivity_specificity_test, 1, min))
  #assigning the identified threshold from min-max to the cutoff variable
  cutoff <- auc_test$thresholds[indx]
  #assign the sensitivity, specificity, and the cutoff for each for each identified model of every alpha value to sn_sp_threshold_matrix_train
  sn_sp_threshold_matrix_test[i,] <- c(sensitivity_specificity_test[indx,][1],sensitivity_specificity_test[indx,][2], cutoff, alphas[i])
  
}

#Using the min-max approach to find the optimal model (alpha value) assosiated with the best balance of sensitivity and specificity on the training and test set
sn_sp_threshold_matrix_train[which.max(apply(sn_sp_threshold_matrix_train, 1, min)),]
sn_sp_threshold_matrix_test[which.max(apply(sn_sp_threshold_matrix_test, 1, min)),]


###Visualize AUC for alpha = 1

cvx_lasso <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = 1, type.measure = "deviance", foldid = foldid)

pred_test_lasso <- predict(cvx_lasso, newx = x[test_index,], type = "response", s=cvx_lasso$lambda.min )[,1]

pred_train_lasso <- predict(cvx_lasso, newx = x[-test_index,], type = "response", s=cvx_lasso$lambda.min )[,1]


auc_train_lasso <- roc(y[-test_index], pred_train_lasso)

auc_test_lasso <- roc(y[test_index], pred_test_lasso)

par(mfrow = c(1,2))

plot(auc_train_lasso)

plot(auc_test_lasso)

par(mfrow = c(1,1))

plot(cvx_lasso)


###coefficients 
coef_lasso <- coef(cvx_lasso, s = cvx_lasso$lambda.min)[,1]
coef_lasso <- coef_lasso[coef_lasso != 0]
coef_lasso_order <- coef_lasso[order(abs(coef_lasso), decreasing = TRUE)]
coef_lasso_order


###looking closer at age
covid_data$Severirty <- as.factor(covid_data$Severirty)

age_severity_logistical <- glm(Severirty ~ AGE, family = "binomial", data = covid_data)
summary(age_severity_logistical)

deviance_residuals <- resid(age_severity_logistical, type = "deviance")

#visualize
table(covid_data$AGE)
table(covid_data$Severirty)
boxplot(AGE ~ Severirty, data = covid_data)

ggplot(covid_data, aes(x = covid_data$Severirty, y = covid_data$AGE)) +
  geom_violin()

#immune status better than age?

###compare to IL-6
covid_data$Severirty <- as.factor(covid_data$Severirty)

il6_severity_logistical <- glm(Severirty ~ covid_data$`IL-6`, family = "binomial", data = covid_data)
summary(il6_severity_logistical)

deviance_residuals <- resid(age_severity_logistical, type = "deviance")
