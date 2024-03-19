###Assignment 3
### Jacob Hambly and Robin Zutshi


library(readxl)
library(ggplot2)
library(GGally)
library(glmnet)
library(pROC)
library(plotly)
library(ggfortify)
library(ggalt)


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
foldid

#train model - 10 folds
alphas <- seq(0,1,0.01)

auc_df <- data.frame(matrix(NA, nrow = length(alphas), ncol = 3))
colnames(auc_df) <- c("Training_AUC", "Testing_AUC", "Alpha_Values")

training_df_10 <- data.frame(alpha = numeric(), log_lambda = numeric(),lambda = numeric(), deviance = numeric(), auc = numeric())

#for each alpha value we get the model corresponding to lambda.min. We pick the lambda.min model as our 'best model for each alpha value. From there we can measure how well each model for each alpha value performs on the training and testing data. So there will be 101 models. For each of the 101 models we also pick the best threshold. The best performing out of the 101 models will be dubbed 'the best model' and have a corresponding threshold.

#10-fold
for (i in 1:length(alphas)) {
  #training the model
  cvx <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = alphas[i], type.measure = "deviance", foldid = foldid)
  
  #storing lambdas, deviance, and alpha values in variables
  lambdas <- cvx$lambda
  deviance <- cvx$cvm
  alpha_values <- rep(alphas[i], length(lambdas))
  auc <- rep(NA, length(lambdas))
  
  #creating a temporary data frame to store the values for alpha, lambda, and deviance
  temp_df <- data.frame(alpha = as.numeric(alpha_values), log_lambda = as.numeric(log(lambdas)), lambda = lambdas,  deviance = as.numeric(deviance), AUC = auc)
  
  #appending temp_df to training_df
  training_df_10 <- rbind(training_df_10, temp_df)
  
  #using the trained model to make predictions of the training data
  pred_train <- predict(cvx, newx = x[-test_index,], type = "response", s=cvx$lambda.min)
  #using the trained model to make predictions for the testing data 
  pred_test <- predict(cvx, newx = x[test_index,], type = "response", s=cvx$lambda.min)
  #calculating the area under the curve for the training predictions 
  auc_train <- roc(y[-test_index], pred_train)
  #calculating the area under the curve for the testing predictions 
  auc_test <- roc(y[test_index], pred_test)
  #appending the area under of the curve for each model to auc_df
  auc_df[i,] <- c(auc_train$auc[1], auc_test$auc[1], alphas[i])
  
}


plot_ly(data =training_df_10, x = ~alpha, y = ~log_lambda, z = ~deviance, color = ~deviance ) %>%
  layout(title = "3D plot of training for 10-fold cross validation")

#training the model with alpha = 0.7 and number of folds = 10
cvx_10fold <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = 0.71, type.measure = "deviance", foldid = foldid) 

#plotting binomial deviance and log(lambda) for alpha = 0.71
plot(cvx_10fold)
  
#predicting using the testing data
pred_test_10fold <- predict(cvx_10fold, newx = x[test_index,], type = "response", s=cvx_10fold$lambda.min)

#predicting using the training data
pred_train_10fold <- predict(cvx_10fold, newx = x[-test_index,], type = "response", s=cvx_10fold$lambda.min )


#finding the AUC for the training data
auc_train_10fold <- roc(y[-test_index], pred_train_10fold)

#finding the AUC for the testing data
auc_test_10fold <- roc(y[test_index], pred_test_10fold)


par(mfrow = c(1,1))



###Looking at the optimal threshold - test

#creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
sn_sp_threshold_test_10fold <- cbind(auc_test_10fold$sensitivities, auc_test_10fold$specificities)

indx_test_10fold <- which.max(apply(sn_sp_threshold_test_10fold, 1, min))
sn_sp_threshold_test_10fold[indx_test_10fold,]

#threshold for the min max sn-sp
threshold_test <- auc_test_10fold$thresholds[indx_test_10fold]


#optimal sn-sp and threshold for 10fold
print(sn_sp_threshold_test_10fold[indx_test_10fold,])
print(threshold_test)


###Looking at the optimal threshold - train 

#creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
sn_sp_threshold_train_10fold <- cbind(auc_train_10fold$sensitivities, auc_train_10fold$specificities)

indx_train_10fold <- which.max(apply(sn_sp_threshold_train_10fold, 1, min))
sn_sp_threshold_train_10fold[indx_train_10fold,]

#threshold for the min max sn-sp
threshold_train <- auc_train_10fold$thresholds[indx_train_10fold]

#optimal sn-sp and threshold for 10fold
sn_sp_threshold_train_10fold[indx_train_10fold,];threshold_train 
sn_sp_threshold_train_10fold[indx_train_10fold,][1]

###add lines for sp and sens and titles for alpha
par(mfrow = c(1,2))

#plotting the AUC plot for the training data
plot(auc_train_10fold)
abline(h = sn_sp_threshold_train_10fold[indx_train_10fold,][1], v = sn_sp_threshold_train_10fold[indx_train_10fold,][2], col = 'blue', lty=2  )

#plotting the AUC curve for the testing data
plot(auc_test_10fold)
abline(h = sn_sp_threshold_test_10fold[indx_test_10fold,][1], v = sn_sp_threshold_test_10fold[indx_test_10fold,][2], col = 'blue', lty=2  )

par(mfrow = c(1,1))


###coefficients - 10 fold

alphas_10_coef <- c(0.7, 0.71)

for (i in 1:length(alphas_10_coef)) {
  
  cvx_10fold <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 10, family = "binomial", alpha = alphas[i], type.measure = "deviance", foldid = foldid)
  
  #extracting the coefficients for alpha=0.45, lambda = lambda.min
  coef_10fold <- coef(cvx_10fold, s = cvx_10fold$lambda.min)[,1]
  #subset all the coefficents  that don't equal to 0
  coef_10fold <- coef_10fold[coef_10fold != 0]
  #order the coefficients in decreasing order by their absolute magnitudes
  coef_10fold_order <- coef_10fold[order(abs(coef_10fold), decreasing = TRUE)]
  print(alphas_10_coef[i])
  print(coef_10fold_order)
}

#training model 20 -folds

#generate fold memberships for 20 fold
set.seed(1717)
cvxf <- cv.glmnet(x[-test_index,],y[-test_index], nfolds = 20, keep = TRUE )
table(cvxf$foldid)
foldid_20 = cvxf$foldid

auc_df_20 <- data.frame(matrix(NA, nrow = length(alphas), ncol = 3))
colnames(auc_df_20) <- c("Training_AUC", "Testing_AUC", "Alpha_Values")

training_df_20 <- data.frame(alpha = numeric(), log_lambda = numeric(),lambda = numeric(), deviance = numeric())


#20 fold
for (i in 1:length(alphas)) {
  #training the model
  cvx <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 20, family = "binomial", alpha = alphas[i], type.measure = "deviance", foldid = foldid_20)
  
  #storing lambdas, deviance, and alpha values in variables
  lambdas <- cvx$lambda
  deviance <- cvx$cvm
  alpha_values <- rep(alphas[i], length(lambdas))
  
  #creating a temporary data frame to store the values for alpha, lambda, and deviance
  temp_df <- data.frame(alpha = as.numeric(alpha_values), log_lambda = as.numeric(log(lambdas)), lambda = lambdas,  deviance = as.numeric(deviance))
  
  #appending temp_df to training_df
  training_df_20 <- rbind(training_df_20, temp_df)
  
  
  #using the trained model to make predictions of the training data
  pred_train_20 <- predict(cvx, newx = x[-test_index,], type = "response", s=cvx$lambda.min)
  
  #using the trained model to make predictions for the testing data 
  pred_test_20 <- predict(cvx, newx = x[test_index,], type = "response", s=cvx$lambda.min )
  
  #calculating the area under the curve for the training predictions 
  auc_train_20 <- roc(y[-test_index], pred_train_20)
  #calculating the area under the curve for the testing predictions 
  auc_test_20 <- roc(y[test_index], pred_test_20)
  #appending the area under of the curve for each model to auc_df
  auc_df_20[i,] <- c(auc_train_20$auc[1], auc_test_20$auc[1], alphas[i])
  
}

plot_ly(data =training_df_20, x = ~alpha, y = ~log_lambda, z = ~deviance, color = ~deviance ) %>%
  layout(title = "3D plot of training for 20-fold cross validation")

#closer look at alpha 0.64
cvx_20fold <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 20, family = "binomial", alpha = 0.64, type.measure = "deviance", foldid = foldid_20)

lambdas_20 <- cvx_20fold$lambda.min

plot(glmnet(x[-test_index,], y[-test_index], alpha = 0.64, lambda = lambdas_20))
  
  #predicting using the testing data
  pred_test_20fold <- predict(cvx_20fold, newx = x[test_index,], type = "response", s=cvx_20fold$lambda.min )[,1]
  
  #predicting using the training data
  pred_train_20fold <- predict(cvx_20fold, newx = x[-test_index,], type = "response", s=cvx_20fold$lambda.min )[,1]
  
  #plotting bionomial deviance and log(lambda) for alpha = 0.64
  plot(cvx_20fold)
  
  
  #finding the AUC for the testing data
  auc_test_20fold <- roc(y[test_index], pred_test_20fold)
  
  #finding the AUC for the training data
  auc_train_20fold <- roc(y[-test_index], pred_train_20fold)
  
  par(mfrow = c(1,1))
  
  
  
  ###Looking at the optimal threshold - test
  
  #creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
  sn_sp_threshold_test_20fold <- cbind(auc_test_20fold$sensitivities, auc_test_20fold$specificities)
  
  indx_test_20fold <- which.max(apply(sn_sp_threshold_test_20fold, 1, min))
  sn_sp_threshold_test_20fold[indx_test_20fold,]
  
  #threshold for the min max sn-sp
  threshold_test_20 <- auc_test_20fold$thresholds[indx_test_20fold]
  

  ###Looking at the optimal threshold - train
  
  #creating a matrix of the sensitivities and specificity for each of the thresholds used to generate the ROC curve
  sn_sp_threshold_train_20fold <- cbind(auc_train_20fold$sensitivities, auc_train_20fold$specificities)
  
  indx_train_20fold <- which.max(apply(sn_sp_threshold_train_20fold, 1, min))
  sn_sp_threshold_train_20fold[indx_train_20fold,]
  
  #threshold for the min max sn-sp
  threshold_train_20 <- auc_train_20fold$thresholds[indx_train_20fold]

  
  
  ###add lines for sp and sens and titles for alpha
  par(mfrow = c(1,2))
  
  #plotting the AUC plot for the training data
  plot(auc_train_20fold)
  abline(h = sn_sp_threshold_train_20fold[indx_train_20fold,][1], v =sn_sp_threshold_train_20fold[indx_train_20fold,][2], col = 'blue', lty = 2 )
  
  #plotting the AUC curve for the testing data
  plot(auc_test_20fold)
  abline(h = sn_sp_threshold_test_20fold[indx_test_20fold,][1], v =sn_sp_threshold_test_20fold[indx_test_20fold,][2], col = 'blue', lty = 2 )
  




###coefficients - 20 fold

alphas_20_coef <- c(0.38, 0.64)

for (i in 1:length(alphas_20_coef)) {
  
  cvx_20fold <- cv.glmnet(x[-test_index,], y[-test_index], nfolds = 20, family = "binomial", alpha = alphas_20_coef[i], type.measure = "deviance", foldid = foldid_20)

#extracting the coefficients for alpha=0.45, lambda = lambda.min
  coef_20fold <- coef(cvx_20fold, s = cvx_20fold$lambda.min)[,1]
#subset all the coefficents  that don't equal to 0
  coef_20fold <- coef_20fold[coef_20fold != 0]
#order the coefficients in decreasing order by their absolute magnitudes
  coef_20fold_order <- coef_20fold[order(abs(coef_20fold), decreasing = TRUE)]
  print(alphas_20_coef[i])
  print(coef_20fold_order)
}



###looking closer at age

#visualize
table(covid_data$AGE)
table(covid_data$Severirty)
boxplot(AGE ~ Severirty, data = covid_data, main = "Boxplot of Age and Severity", xlab = "Severity", ylab = "Age")


hist(covid_data$AGE, xlab = "Age", main = "Histogram of age of patients")

covid_data$Severirty <- ifelse(covid_data$Severirty == "Mild", 0, 1)

age_severity_logistical <- glm(Severirty ~ AGE, family = "binomial", data = covid_data)

summary(age_severity_logistical)

coef(age_severity_logistical) 

plot(covid_data$AGE, covid_data$Severirty)

#Generating table of age and severity 
table_prop_severe <- table(covid_data$AGE, covid_data$Severirty)

proportion_severe <- table_prop_severe[,2]/apply(table_prop_severe, 1, sum)

proportion_severe2 <- ifelse(proportion_severe==1, .99, proportion_severe)

proportion_severe3 <- ifelse(proportion_severe2==0, .01, proportion_severe2)

odds_severe <- proportion_severe3/(1-proportion_severe3)

log_odds_severe <- log(odds_severe)  

plot(as.numeric(dimnames(table_prop_severe)[[1]]), log_odds_severe, xlab = "Age", ylab = "Log odds of Severe COVID-19", main = "Log odds of Severe COVID-19 and Age")
abline(glm(Severirty ~ AGE, family = "binomial", data = covid_data))



#Point-Biserial Correlation
cor_age_sev <- cor.test(covid_data$AGE, as.numeric(covid_data$Severirty))


##### Problem 2


### Data acquisition and formatting

# Load objects saved in .rda file into global environment
load(file = "geneexpression2.rda")

# Rename new data frame into something easier to track
gene_expression_data <- dat

# Remove original object loaded from file
rm(dat)

# View new, renamed data frame
View(gene_expression_data)


### Data exploration

# Check the dimensions o the new data frame
dim(gene_expression_data)

# Check out the structure of the data frame
str(gene_expression_data)

# Quick summary statistics of dataset
summary(gene_expression_data)

#seeing if there are any duplicate rows
dim(gene_expression_data);dim(unique(gene_expression_data))

#seeing if there are any duplicate columns
dim(gene_expression_data);dim(unique(gene_expression_data, MARGIN = 2))


# Check the mean of data 
apply(gene_expression_data, 2, mean)

# Check the standard deviation of data
apply(gene_expression_data, 2, sd)

# Extracting the first four columns of the original data frame and combining them with the standardized values of the remaining columns
gene_expression_standardized <- as.data.frame(scale(gene_expression_data))

# Check to see if mean of data is 0
apply(gene_expression_standardized, 2, mean)

# Check to see if standard deviation of data is 1
apply(gene_expression_standardized, 2, sd)

# View standardized data frame
View(gene_expression_standardized)


### Metadata extraction

# Add a column to data frame to specify if cell types are Naive (NAI), Effector (EFFE), or Memory (MEM) cells
gene_expression_standardized$cell_type <- NA

# Identify and label 'Naive' cells based on row names containing "NAI"
gene_expression_standardized$cell_type[grepl("NAI", rownames(gene_expression_standardized))] <- "Naive"

# Identify and label 'Effector' cells based on row names containing "EFFE"
gene_expression_standardized$cell_type[grepl("EFFE", rownames(gene_expression_standardized))] <- "Effector"

# Identify and label 'Memory' cells based on row names containing "MEM"
gene_expression_standardized$cell_type[grepl("MEM", rownames(gene_expression_standardized))] <- "Memory"

# Change class of 'cell_type' column to factor 
gene_expression_standardized$cell_type <- as.factor(gene_expression_standardized$cell_type)

# Add a column to data frame to specify if observations are of Healthy (HEA) or Melanoma (MEL) status
gene_expression_standardized$status <- NA

# Identify and label observations as 'Healthy' based on row names containing "HEA"
gene_expression_standardized$status[grepl("HEA", rownames(gene_expression_standardized))] <- "Healthy"

# Identify and label observations as 'Melanoma' based on row names containing "MEL"
gene_expression_standardized$status[grepl("MEL", rownames(gene_expression_standardized))] <- "Melanoma"

# Change class of 'status' column to factor 
gene_expression_standardized$status <- as.factor(gene_expression_standardized$status)

# Rearrange columns such that the last two columns are now the first two
gene_expression_standardized <- gene_expression_standardized[c((ncol(gene_expression_standardized)-1):ncol(gene_expression_standardized), 1:(ncol(gene_expression_standardized)-2))]

# Check if changes happened in data frame
View(gene_expression_standardized)

# Check class of 'cell_type' and 'status' columns
class(gene_expression_standardized$cell_type)
class(gene_expression_standardized$status)


### Data analysis

# Correlation plot analysis of numerical variables
#ggpairs(gene_expression_data[,c(-1,-2)])

# Create a covariance matrix of the numeric variables
covar_matrix <- var(gene_expression_standardized[,c(-1,-2)])

# Round covariance matrix to 3 sf.
round(covar_matrix, 3)

# Calculate eigenvalues and eigenvectors of data covariance matrix
eigen <- eigen(covar_matrix)

# Store calculated eigenvalues of data matrix into new list
eigenvalues <- eigen$values

# Store calculated eigenvectors of data matrix into new object
eigenvectors <- eigen$vectors

# Object to store number of columns needed to create result matrix
ncol = length(eigenvectors[1,])

# Object to store number of rows needed to create result matrix
nrow = (length(eigenvectors[1,]) + 2)

# Create a matrix to hold results, eigenvalues, and cumulative variances
result_matrix <- matrix(NA, ncol = ncol, nrow = nrow)

# Assign eigenvalues to result_matrix
result_matrix[nrow(result_matrix) - 1, ] <- eigenvalues

# Cumulative sum of PC variances 
result_matrix[nrow(result_matrix), ] <- cumsum(eigenvalues/sum(eigenvalues)*100)

# Assign row names to result_matrix: gene names, eigenvalue estimate, cumulative percentage
rownames(result_matrix) <- c(colnames(gene_expression_standardized[,c(-1,-2)]), 
                             "lambda.hat", "Cumulative %")

# Assign column names to result_matrix: PC1, PC2, ..., PCn
colnames(result_matrix) <- paste("PC", 1:ncol(result_matrix), sep = "")

# Fill the top rows of result_matrix with eigenvectors
result_matrix[1:(nrow(result_matrix) - 2), 1:ncol(result_matrix)] <- eigenvectors

# Round result matrix to 3 sf.
round(result_matrix, 3)

# Having a look at the first two PCs
result_matrix[1:(nrow(result_matrix) - 2), 1:2]


### Visualization plots

# Scree plot of PCs vs their eigenvalues
plot(x = 1:length(eigenvalues), y = eigenvalues, type = "b", 
     ylab = "Eigenvalue", xlab = "Component", main = "Scree Plot")

# Plot of PCs vs cumulative variance explained
plot(x = 1:length(eigenvalues), y = result_matrix[nrow(result_matrix),], type = "b", 
     ylab = "Cumulative Variance (%)", xlab = "Component", main = "Cumulative Variance Plot")


# Create an empty matrix with as many rows as observations and as many cols as number of variables in original dataset
pc_scores <- matrix(NA, nrow = nrow(gene_expression_standardized), 
                    ncol = ncol(gene_expression_standardized[,-c(1:2)]))

# View new matrix to hold PC scores
pc_scores

# Write a function to calculate the PC score for given eigenvector
score_fn <- function (i, eigenvector) { # Here 'i' is the 'i-th' observation in original dataset
  # Perform matrix multiplication between transpose of given eigenvector and observation 'i'  
  as.numeric(t(eigenvector) %*% i)
}

# Compute scores for all PCs
for (j in 1:ncol(eigenvectors)) {  # Here j is the j-th PC
  # Apply the score function to compute the scores for the j-th principal component
  pc_scores[, j] <- apply(gene_expression_standardized[,-c(1:2)], MARGIN = 1, score_fn, 
                          eigenvector = eigenvectors[, j]) 
}

# Check covariance between PC scores
var(pc_scores)

# Compare diagonal variance of PC scores and eigenvalues
all.equal(diag(var(pc_scores)), eigenvalues)

# Computation shows diagonal covariance of scores is equal to eigenvalues


# Generate a PCA plot using calculated PC scores
ggplot(data = pc_scores) +
  geom_point(mapping = aes(x = pc_scores[, 1], y = pc_scores[, 2], 
                           color = gene_expression_standardized$cell_type,
                           shape = gene_expression_standardized$status)) +
  labs(title = "PCA Plot", 
       x = paste("PC1"," ","(",round(result_matrix[nrow(result_matrix), 1], 2),"%)", sep = ""), 
       y = paste("PC2"," ","(",round(result_matrix[nrow(result_matrix), 2] - 
                                       result_matrix[nrow(result_matrix), 1], 2),"%)", sep = ""),
       color = "Cell Type",
       shape = "Status") +
  geom_encircle(aes(x = pc_scores[, 1], y = pc_scores[, 2], 
                    group = gene_expression_standardized$cell_type), 
                expand = 0.02, color = "black", alpha = 0.3)



# Create a biplot using autoplot and prcomp with colored PC scores by cell type
autoplot(prcomp(gene_expression_standardized[,-c(1:2)]), data = gene_expression_standardized, colour = 'cell_type', shape = 'status', loadings = TRUE, loading.lables = TRUE, loadings.colour = 'darkgrey') +
  labs(title = "PCA Plot")


