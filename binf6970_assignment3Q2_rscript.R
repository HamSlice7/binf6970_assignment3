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
