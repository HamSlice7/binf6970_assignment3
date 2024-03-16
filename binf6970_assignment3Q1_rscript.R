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
x$SEX <- ifelse(x$SEX == "M", 1, 0)
x <- as.matrix(x)

y <- ifelse(covid_data_standardized$Severirty == "Mild", 0, 1)
table(y)

###Train-test
set.seed(1717)

#generating indexes for the test set
test_index <- sample(dim(x)[1],round(dim(x)[1] *.25), replace = FALSE )

#generate fold memberships for 10 fold
cvxf <- cv.glmnet(x[-test_index,],y[-test_index], nfolds = 10, keep = TRUE )
table(cvxf$foldid)
foldid = cvxf$foldid

#train model - 10 folds
alphas <- seq(0,1,0.01)

auc_df <- data.frame(matrix(NA, nrow = length(alphas), ncol = 3))
colnames(auc_df) <- c("Training_AUC", "Testing_AUC", "Alpha_Values")

#for each alpha value we get the model corresponding to lambda.min. We pick the lambda.min model as our 'best model for each alpha value. From there we can measure how well each model for each alpha value performs on the training and testing data. So there will be 101 models. For each of the 101 models we also pick the best threshold. The best performing out of the 101 models will be dubbed 'the best model' and have a corresponding threshold.

#10-fold
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
  
}

#getting the model associated with the largest AUC
auc_df[which.max(auc_df$Testing_AUC),]

##taking a closer look at the model with the best AUC

#training the model with alpha = 0.7 and number of folds = 10
cvx_10fold <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = 0.7, type.measure = "deviance", foldid = foldid)

#predicting using the testing data
pred_test_10fold <- predict(cvx_10fold, newx = x[test_index,], type = "response", s=cvx_10fold$lambda.min )[,1]

#predicting using the training data
pred_train_10fold <- predict(cvx_10fold, newx = x[-test_index,], type = "response", s=cvx_10fold$lambda.min )[,1]

#plotting bionomial deviance and log(lambda) for alpha = 0.7
plot(cvx_10fold)

#finding the AUC for the training data
auc_train_10fold <- roc(y[-test_index], pred_train_10fold)

#finding the AUC for the testing data
auc_test_10fold <- roc(y[test_index], pred_test_10fold)

par(mfrow = c(1,2))

#ploting the AUC plot for the training data
plot(auc_train_10fold)

#plotting the AUC curve for the testing data
plot(auc_test_10fold)

par(mfrow = c(1,1))



###Looking at the optimal threshold - test

#creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
sn_sp_threshold_test_10fold <- cbind(auc_test_10fold$sensitivities, auc_test_10fold$specificities)

indx_test_10fold <- which.max(apply(sn_sp_threshold_test_10fold, 1, min))
sn_sp_threshold_test_10fold[indx_test_10fold,]

#threshold for the min max sn-sp
threshold_test <- auc_test_10fold$thresholds[indx_test_10fold]


#optimal sn-sp and threshold for 10fold
sn_sp_threshold_test_10fold[indx_test_10fold,];threshold_test 



###Looking at the optimal threshold - train

#creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
sn_sp_threshold_train_10fold <- cbind(auc_train_10fold$sensitivities, auc_train_10fold$specificities)

indx_train_10fold <- which.max(apply(sn_sp_threshold_train_10fold, 1, min))
sn_sp_threshold_train_10fold[indx_train_10fold,]

#threshold for the min max sn-sp
threshold_train <- auc_train_10fold$thresholds[indx_train_10fold]

#optimal sn-sp and threshold for 10fold
sn_sp_threshold_train_10fold[indx_train_10fold,];threshold_train 


###coefficients - 10 fold
#extracting the coefficients for alpha=0.7, lambda = lambda.min
coef_10fold <- coef(cvx_10fold, s = cvx_10fold$lambda.min)[,1]
#subset all the coefficents  that don't equal to 0
coef_10fold <- coef_10fold[coef_10fold != 0]
#order the coefficients in decreasing order by their absolute magnitudes
coef_10fold_order <- coef_10fold[order(abs(coef_10fold), decreasing = TRUE)]
coef_10fold_order

#training model 20 -folds

#generate fold memberships for 20 fold
cvxf <- cv.glmnet(x[-test_index,],y[-test_index], nfolds = 20, keep = TRUE )
table(cvxf$foldid)
foldid_20 = cvxf$foldid

auc_df_20 <- data.frame(matrix(NA, nrow = length(alphas), ncol = 3))
colnames(auc_df_20) <- c("Training_AUC", "Testing_AUC", "Alpha_Values")

#20 fold
for (i in 1:length(alphas)) {
  #training the model
  cvx <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 20, family = "binomial", alpha = alphas[i], type.measure = "deviance", foldid = foldid_20)
  #using the trained model to make predictions of the training data
  pred_train_20 <- predict(cvx, newx = x[-test_index,], type = "response", s=cvx$lambda.min)[,1]
  #using the trained model to make predictions for the testing data 
  pred_test_20 <- predict(cvx, newx = x[test_index,], type = "response", s=cvx$lambda.min )[,1]
  #calculating the area under the curve for the training predictions 
  auc_train_20 <- roc(y[-test_index], pred_train_20)
  #calculating the area under the curve for the testing predictions 
  auc_test_20 <- roc(y[test_index], pred_test_20)
  #appending the area under of the curve for each model to auc_df
  auc_df_20[i,] <- c(auc_train_20$auc[1], auc_test_20$auc[1], alphas[i])
  
}

#getting the model associated with the largest AUC on the testing data
auc_df_20[which.max(auc_df_20$Testing_AUC),]

##taking a closer look at the model with the best AUC

#training the model with alpha = 0.45 and number of folds = 20
cvx_20fold <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 20, family = "binomial", alpha = 0.45, type.measure = "deviance", foldid = foldid_20)

#plotting bionomial deviance and log(lambda) for alpha = 0.45
par(mfrow = c(1,1))
plot(cvx_20fold)

#predicting using the testing data
pred_test_20fold <- predict(cvx_20fold, newx = x[test_index,], type = "response", s=cvx_20fold$lambda.min )[,1]

#predicting using the training data
pred_train_20fold <- predict(cvx_20fold, newx = x[-test_index,], type = "response", s=cvx_20fold$lambda.min )[,1]


#finding the AUC for the training data
auc_train_20fold <- roc(y[-test_index], pred_train_20fold)

#finding the AUC for the testing data
auc_test_20fold <- roc(y[test_index], pred_test_20fold)

par(mfrow = c(1,2))

#ploting the AUC plot for the training data
plot(auc_train_20fold)

#plotting the AUC curve for the testing data
plot(auc_test_20fold)

par(mfrow = c(1,1))

###Looking at the optimal threshold - test

#creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
sn_sp_threshold_test_20fold <- cbind(auc_test_20fold$sensitivities, auc_test_20fold$specificities)

indx_test_20fold <- which.max(apply(sn_sp_threshold_test_20fold, 1, min))
sn_sp_threshold_test_20fold[indx_test_20fold,]

#threshold for the min max sn-sp
threshold_test_20fold <- auc_test_20fold$thresholds[indx_test_20fold]


#optimal sn-sp and threshold for 20 fold
sn_sp_threshold_test_20fold[indx_test_20fold,];threshold_test_20fold


###Looking at the optimal threshold - training

#creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
sn_sp_threshold_train_20fold <- cbind(auc_train_20fold$sensitivities, auc_train_20fold$specificities)

indx_train_20fold <- which.max(apply(sn_sp_threshold_train_20fold, 1, min))
sn_sp_threshold_train_20fold[indx_train_20fold,]

#threshold for the min max sn-sp
threshold_train_20fold <- auc_train_20fold$thresholds[indx_train_20fold]


#optimal sn-sp and threshold for 20 fold
sn_sp_threshold_train_20fold[indx_train_20fold,];threshold_train_20fold


###coefficients - 10 fold
#extracting the coefficients for alpha=0.45, lambda = lambda.min
coef_20fold <- coef(cvx_20fold, s = cvx_20fold$lambda.min)[,1]
#subset all the coefficents  that don't equal to 0
coef_20fold <- coef_20fold[coef_20fold != 0]
#order the coefficients in decreasing order by their absolute magnitudes
coef_20fold_order <- coef_20fold[order(abs(coef_20fold), decreasing = TRUE)]
coef_20fold_order


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
