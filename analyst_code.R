# Part 4: Analysis
################# Noise filtering & Dimentionality reduction ###############
## load libraries
library(cluster)
library(dplyr)
library(tidyr)
library(tibble)
library(gplots)

## import data
data<-read.csv("expression_data.csv", row.names=1)

# 1) Expressed in at least 20% of samples (i.e. for each gene, at least 20% of the gene-expression values must be > log2(15)
# result: 39661 genes passed 
data1 <- rowSums(data > log2(15)) >= 0.2*ncol(data)
exp_filter<- data[data1, ]
# check how many genes are left 
print(c("Number of probes left after filter 1:",nrow(exp_filter)))

# 2) Have a variance significantly different from the median variance of all probe sets using a threshold of p<0.01
# result: 15508 genes passed chisq filter

## calculate row variance
rvar<- apply(exp_filter, 1, var)
#exp_filter$Variance <- rvar
mvar<- median(rvar)
df<- ncol(data) - 1

## filter: genes that pass the one-tail chi sq test 
chi_filter <- exp_filter[133*(rvar/mvar)>qchisq(0.99,133),]
## check how many genes left
print(c("Number of probes left after filter2:",nrow(chi_filter)))

# 3) Have a coefficient of variation > 0.186.
# result: 1531 genes passed var_filter
cv<-apply(chi_filter, 1, sd)/ apply(chi_filter, 1, mean)
var_filter<-chi_filter[cv> 0.186,]
print(c("Number of probes left after filter3:", nrow(var_filter)))

# 4) Write out a different file containing the gene expression matrix for genes passing all three of the filters from 4.1, 4.2, and 4.3.
write.csv(var_filter, sep='', file="data4.4.csv")

# 5) For groups with Biologist role only: Write out the expression matrix for probesets that pass the expression threshold from 4.2 to a file with write.csv.
write.csv(chi_filter,sep = "",file = "data4.5.csv")


################# Part 5: Hierarchical clustering & subtype discovery ##################
# 1) Perform hierarchical clustering on your fully filtered data matrix from Part 4.4. Be sure to check that you are clustering the patients and not the genes.
# use hclust() and euclidean distance method for clustering
clusters<- hclust(dist(t(var_filter)))
# plot dendrogram
plot(clusters,labels = FALSE,main = "Cluster Dendrogram",sub="")

# 2) Cut the dendrogram such that the samples are divided into 2 clusters. How many samples are in each cluster?
## cut dendrogram into 2 clusters
#result: 134 genes in cluster1, 2 genes in cluster2
cut_avg<- cutree(clusters, k=2)
table(clusters)
print(c('Number of samples in cluster 1: ', sum(cut_avg==1)))
print(c('Number of samples in cluster 2: ', sum(cut_avg==2)))

# 3) Create a heatmap of the gene-expression of each gene across all samples using the heatmap() function. Include a column colorbar by setting the ColSideColors variable in the heatmap function equal to "red" if the sample belongs to the C3 subtype and "blue" otherwise. Subtype annotation can be found in the annotation matrix under the title cit-coloncancermolecularsubtype.

## load metadata
proj_metadata<- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
metadata<- subset(proj_metadata, select=c(geo_accession, cit.coloncancermolecularsubtype))
target<- metadata$geo_accession
# go through data vector and for each matching C3 subtype, append blue
C3_subtype<- c()
  for (i in 1:col(metadata)){
    C3_subtype[i] <- substr(colnames(metadata)[i],1,9)
}
C3_subtype<- as.factor(C3_subtype)
C3_subtype == target
## create heatmap 
color<-ifelse(proj_metadata$cit.coloncancermolecularsubtype == 'C3','red','blue')
heatmap(as.matrix(var_filter),ColSideColors = color, main = '', xlab='Samples',
        ylab='Genes')

# 4) Using the expression matrix from Part 4.4 and the cluster memberships from Part 5.2, identify genes differentially expressed between the two clusters using a Welch t-test (results in a ranked list). Write out a dataframe containing the probeset ID, t-statistic, p-value, and adjusted p-value (i.e. FDR, see the p.adjust function) to a comma separated file for each comparison. How many genes are differentially expressed at adjusted p<0.05 between the clusters for both lists?
## separate data into cluster1 and cluster2 for t-test  
# cluster 1
group1<-var_filter[, cut_avg==1]
# cluster 2
group2<-var_filter[, cut_avg==2]
## Welch t-test (does now assume that the two groups have equal variance)
## alternative hypothesis: true difference in means is not equal to 0

## t.test(group1, group2, alternative=c("two-sided"))
t_data<-apply(as.matrix(var_filter),1,function(x) t.test(x=x[cut_avg==1],y=x[cut_avg==2]))
## extract p-val and test statistics
p_value<- sapply(t_data,function(x) x$p.value)
test_stat<- sapply(t_data,function(x) x$statistic)
## compute q-value 
adj_p<- p.adjust(p_value,method = "fdr")
## build new dataframe with probeset ID, test statistic, p-val and q-val
data5.4<- data.frame("Probeset_ID" = c(row.names(var_filter), test_stat,p_value,adj_p))
print(data5.4)

# 5) Select the most differentially expressed genes that you feel best define the clusters and explain your selection.
# number of significant (p < 0.5) genes is 1249 
sig_genes<- filter(data5.4, p_value < 0.05)
print(c("The number of significant (p < 0.5) genes is: ", nrow(sig_genes)))
## select those with q-val less than 0.05
differential<- data5.4$Probeset_ID[data5.4$adj_p<0.05]
print(c("Number of probes left with q-val<0.05:", length(differential)))
## save data into csv file
write.csv(data5.4,sep = ",", row.names = FALSE,file = "data5.4.csv")

# 6) For groups with Biologist role only: perform the t-test analysis described in 5.4 on the expression matrix from 4.5 as well and provide to the Biologist as a csv file.
## divide into 2 clusters again
group1_2<-chi_filter[, cut_avg==1]
group2_2<-chi_filter[, cut_avg==2]
# implement t-test analysis on expression matrix from 4.5
for (i in 1:nrow(group1_2)) {
  ttest<- t.test(group1_2[i,], group2_2[i,])
  data5.4[i,]$test_stat <- ttest$statistic
  data5.4[i,]$p_val <- ttest$p.value
}
data5.4$adjust_p <- p.adjust(data5.4$p_val, method = "fdr")
write.csv(data5.4, "data5.6.csv")
