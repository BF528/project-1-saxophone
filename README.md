# Microarray Based Tumor Classification

The purpose of this project is to reproduce analysis of microarray data comparing C3 and C4 colon cancer tumor subtypes from Marisa et. al. Data obtained from GEO Accession #GSE39582

Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

# Contributors

Sooyoun Lee leesu@bu.edu - Programmer

Daniel Goldstein djgoldst@bu.edu - Data Curator

Jason Rose jjrose@bu.edu - Biologist

Sunny Yang yang98@bu.edu - Analyst

# Repository Contents

Programmer.R - Data pre-processing and quality control script using the RMA algorithm to normalize all microarrays together, ComBat to correct for batch effects, and PCA to visualize the quality of the sample distribution.

analyst_code.R - Data processing using noise filtering and dimension reduction, and data analysis using unsupervised hierarchical clustering to identify differential gene expression between the two clusters.
