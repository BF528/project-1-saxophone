# install packages 'Bioconductor' and some libraries(R 4.0.3 version)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install(c("affy","affyPLM","sva",'AnnotationDbi','hgu133plus2.db'))
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)


## 2. Read CEL files and Normalize
# read all 134 gz files
files = list.files("/projectnb/bf528/users/saxophone/project_1/samples/CEL_files", full.names = TRUE)
length(files) # check 134 files

data_affy <- ReadAffy(filenames=files) # read gz files, not unzipping.
data_norm <- rma(data_affy) # rma algorithm

## 3. Plot RLE and NUSE
data_plm <- fitPLM(data_affy, normalize = T, background = T)

result_rle <- RLE(data_plm,type="stats")
median_rle <- result_rle[1,] # extract median results
hist(median_rle, main = 'Histogram of median for rle',
     xlab = 'median for rle',col = 'black', nclass = 10)

result_nuse <- NUSE(data_plm,type="stats")
median_nuse <- result_nuse[1,]
hist(median_nuse, main="Histogram of median RLE", 
     xlab="median RLE", col="black", nclass=10)

## 4. csv file for RMA normalized
norm_exprs <- exprs(data_norm)
write.csv(norm_exprs, 'normdata.csv')

## 5. Batch effects Correction
metadata <- read.csv("/projectnb/bf528/users/saxophone/project_1/samples/proj_metadata.csv")
mod <- model.matrix(~normalizationcombatmod, data = metadata)
combat_adj <- ComBat(dat = norm_exprs, batch = metadata$normalizationcombatbatch, mod = mod)
write.csv(combat_adj,'combat_adj.csv')

## 6. PCA : 2 PCs explain about 20% for total variablity.
combat_adj_t <-t(combat_adj)
combat_scaled <-t(scale(combat_adj_t, center=TRUE,scale=TRUE))
combat_pca <-prcomp(combat_scaled,scale=FALSE,center = FALSE)
summary(combat_pca)
combat_pca$rotation

pca_value <- as.data.frame(combat_pca$rotation)
ggplot(data = pca_value, mapping = aes(x = PC1, y = PC2)) +
  geom_point()+
  theme_bw() +
  labs(title = 'PCA plot', x= 'PC1 : 11.47%', y='PC2 : 8.41%')

## 7. 3 Outliers are in PC2
ggplot(data = pca_value, mapping = aes(y = PC1)) +
  geom_boxplot()+
  theme_bw() +
  labs(title = 'PCA1 Boxplot for detecting outliers')

ggplot(data = pca_value, mapping = aes(y = PC2)) +
  geom_boxplot()+
  theme_bw() +
  labs(title = 'PCA2 Boxplot for detecting outliers')