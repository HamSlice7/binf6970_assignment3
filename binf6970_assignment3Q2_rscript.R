load(file = 'geneexpression2.rda')

head(dat)

dim(dat)

#is the data standardized?

apply(dat, 2, mean)

apply(dat, 2, sd)

#standardizing the data
dat_standardized <- scale(dat)

apply(dat_standardized, 2, mean)

apply(dat_standardized, 2, sd)

#perform pca
pca_ge <- prcomp(dat_standardized)

pca_ge_scores <- pca_ge$x

ggplot(data = pca_ge_scores, aes(x = PC1, y = PC2, col = rownames(pca_ge_scores))) + 
  geom_point()

################################################################################

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


### Metadata extraction

# Add a column to data frame to specify if cell types are Naive (NAI), Effector (EFFE), or Memory (MEM) cells
gene_expression_data$cell_type <- NA

# Identify and label 'Naive' cells based on row names containing "NAI"
gene_expression_data$cell_type[grepl("NAI", rownames(gene_expression_data))] <- "Naive"

# Identify and label 'Effector' cells based on row names containing "EFFE"
gene_expression_data$cell_type[grepl("EFFE", rownames(gene_expression_data))] <- "Effector"

# Identify and label 'Memory' cells based on row names containing "MEM"
gene_expression_data$cell_type[grepl("MEM", rownames(gene_expression_data))] <- "Memory"

# Change class of 'cell_type' column to factor 
gene_expression_data$cell_type <- as.factor(gene_expression_data$cell_type)

# Add a column to data frame to specify if observations are of Healthy (HEA) or Melanoma (MEL) status
gene_expression_data$status <- NA

# Identify and label observations as 'Healthy' based on row names containing "HEA"
gene_expression_data$status[grepl("HEA", rownames(gene_expression_data))] <- "Healthy"

# Identify and label observations as 'Melanoma' based on row names containing "MEL"
gene_expression_data$status[grepl("MEL", rownames(gene_expression_data))] <- "Melanoma"

# Change class of 'status' column to factor 
gene_expression_data$status <- as.factor(gene_expression_data$status)

# Rearrange columns such that the last two columns are now the first two
gene_expression_data <- gene_expression_data[c((ncol(gene_expression_data)-1):ncol(gene_expression_data), 
                                               1:(ncol(gene_expression_data)-2))]

# Check if changes happened in data frame
View(gene_expression_data)

# Check class of 'cell_type' and 'status' columns
class(gene_expression_data$cell_type)
class(gene_expression_data$status)


### Data analysis

# Correlation plot analysis of numerical variables
#ggpairs(gene_expression_data[,c(-1,-2)])

# Create a covariance matrix of the numeric variables
covar_matrix <- var(gene_expression_data[,c(-1,-2)])

# Round covariance matrix to 3 sf.
round(covar_matrix, 3)

# Calculate eigenvalues and eigenvectors of data covariance matrix
eigen <- eigen(covar_matrix)

# Store calculated eigenvalues of data matrix into new list
eigenvalues <- eigen$values

# Store calculated eigenvectors of data matrix into new object
eigenvectors <- eigen$vectors

# 
ncol = length(eigenvectors[1,])

#
nrow = (length(eigenvectors[1,]) + 2)

# Create a matrix to hold results, eigenvalues, and cumulative variances
result_matrix <- matrix(NA, ncol = ncol, nrow = nrow)

# Assign eigenvalues to result_matrix
result_matrix[nrow(result_matrix) - 1, ] <- eigenvalues

# Cumulative sum of PC variances 
result_matrix[nrow(result_matrix), ] <- cumsum(eigenvalues/sum(eigenvalues)*100)

#
rownames(result_matrix) <- c(colnames(gene_expression_data[,c(-1,-2)]), "lambda.hat", "Cumulative %")

#
colnames(result_matrix) <- paste("PC", 1:ncol(result_matrix), sep = "")

#
result_matrix[1:(nrow(result_matrix) - 2), 1:ncol(result_matrix)] <- eigenvectors

# Round result matrix to 3 sf.
round(result_matrix, 3)

# Having a look at the first two PCs
result_matrix[1:(nrow(result_matrix) - 2), 1:2]


### Visualization plots

# Scree plot of PCs vs their eigenvalues
plot(x = 1:length(eigenvalues), y = eigenvalues, type = "b", 
     ylab = "lambda-hat", xlab = "Component", main = " The Scree Plot")

# Plot of PCs vs cumulative variance explained
plot(x = 1:length(eigenvalues), y = result_matrix[nrow(result_matrix),], type = "b", 
     ylab = "Cumulative Variance (%)", xlab = "Component", main = "Cumulative Variance Plot")



#pc_data <- princomp(gene_expression_data[,-c(1:2)], cor = FALSE, scores = TRUE)

#attributes(pc_data)

#pc_scores <- pc_data$scores



# Create an empty matrix with as many rows as observations and as many cols as number of variables in original dataset
pc_scores <- matrix(NA, nrow = nrow(gene_expression_data), ncol = ncol(gene_expression_data[,-c(1:2)]))

# View new matrix to hold PC scores
pc_scores

# Write a function to calculate the PC score for given eigenvector
score_fn <- function (i, eigenvector) { # Here 'i' is the 'i-th' observation in the original dataset
  # Perform matrix multiplication between transpose of given eigenvector and observation 'i'  
  as.numeric(t(eigenvector) %*% i)
}

# Compute scores for all PCs
for (j in 1:ncol(eigenvectors)) {  # Here j is the j-th PC
  # Apply the score function to compute the scores for the j-th principal component
  pc_scores[, j] <- apply(gene_expression_data[,-c(1:2)], MARGIN = 1, score_fn, 
                          eigenvector = eigenvectors[, j]) 
}

# Check covariance between PC scores
var(pc_scores)

# Compare diagonal variance of PC scores and eigenvalues
all.equal(diag(var(pc_scores)), eigenvalues)
# Computation shows diagonal covariance of scores is equal to eigenvalues


# Generate a PCA plot using PC scores
plot(pc_scores[, 1], pc_scores[, 2], main = "PCA Plot", 
     xlab = paste("PC1"," ","(",round(result_matrix[nrow(result_matrix), 1], 2),"%)", sep = ""), 
     ylab = paste("PC2"," ","(",round(result_matrix[nrow(result_matrix), 2] - 
                                        result_matrix[nrow(result_matrix), 1], 2),"%)", sep = ""))

ggplot(data = pc_scores) +
  geom_point(mapping = aes(x = pc_scores[, 1], y = pc_scores[, 2], 
                           color = gene_expression_data$cell_type,
                           shape = gene_expression_data$status)) +
  labs(title = "PCA Plot", 
       x = paste("PC1"," ","(",round(result_matrix[nrow(result_matrix), 1], 2),"%)", sep = ""), 
       y = paste("PC2"," ","(",round(result_matrix[nrow(result_matrix), 2] - 
                                       result_matrix[nrow(result_matrix), 1], 2),"%)", sep = ""),
       color = "Cell Type",
       shape = "Status") +
  geom_encircle(aes(x = pc_scores[, 1], y = pc_scores[, 2], 
                    group = gene_expression_data$cell_type), 
                expand = 0.02, color = "black", alpha = 0.3)

