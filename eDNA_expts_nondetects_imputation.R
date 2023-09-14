setwd("/Users/aholmes/desktop/DSM_eDNA/expts")

library(HTqPCR)
library(mvtnorm)
library(nondetects)

#R version 4.2.1 (2022-06-23)

#read in table with Ct values for 432 technical replicates
#columns are samples (i.e., technical reps of biological reps (individual filters) of an experimental treatment)
#rows are features (in gene expression data features are genes, but here features are used for 6 experimental treatments)
#needs to be a matrix
#code adapted from https://support.bioconductor.org/p/118400/

in_data = read.table("readCt_DSM_expts.txt", 
                     sep = "\t", 
                     stringsAsFactors = FALSE,
                     header = TRUE, 
                     row.names = 1, 
                     na.strings = "NA")

mat = as.matrix(in_data)

#nondetects will need to know which values to impute
#first, create a df with featureCategory with "OK" for all Ct values
feat_cat = as.data.frame(array("OK", dim=dim(mat)), stringsAsFactors = FALSE)

#then, change the featureCategory to "Undetermined" for "NA" entries in the data set
#undetermined values will be imputed later on
feat_cat[is.na(mat)] = "Undetermined" 
rownames(feat_cat) = rownames(mat)
colnames(feat_cat) = colnames(mat)

#boxplot and qPCRimput won't work with NA values so I reimported a version of in_data with NA values set to the Ct max (50)
#there's probably a better way to do this by setting feat_cat values at and above LOD (40.38) to "Undetermined"
in_data = read.table("readCt_DSM_expts_sub50.txt", 
                     sep = "\t", 
                     stringsAsFactors = FALSE,
                     header = TRUE, 
                     row.names = 1, 
                     na.strings = "NA")

mat = as.matrix(in_data)

# Create new instance of qPCRset; needs featureNames to be added
raw = new("qPCRset", exprs = mat, featureCategory = feat_cat)
featureNames(raw) = rownames(mat)

#read in a file with the metadata (filter type and biological replicate) for each sample (i.e., technical replicate)
phenoData <- read.csv("phenoData_14Feb23.csv", header = TRUE)

#add columns for the sample metadata
pData(raw)$sampleName <- phenoData$sampleName
pData(raw)$sampleType <- phenoData$sampleType
pData(raw)$filterType <- phenoData$filterType
pData(raw)$bioRep <- phenoData$bioRep

#set the experimental treatments as target genes
featureNames(raw) <- c("5NTU.500","50NTU.500","PF.500","5NTU.1000","50NTU.1000","PF.1000")
featureType(raw) <- c("Target","Target","Target","Target","Target","Target")

show(raw) #view the qPRCset object named raw
pData(raw) #retrieves sample metadata (i.e., information on experimental phenotypes)
exprs(raw) #retrieves CT values
#each "feature" is one of 6 experimental treatments (turbidity, prefiltration status, and amount of fish tank water added to each bottle)
#the "samples" are technical replicates (individual qPCR wells)

#examine the residuals in the raw qPCR data with max cycles (Ct  50) substituted for samples with undetermined Ct
conds <- paste(pData(raw)$sampleType)

resids <- matrix(nrow=nrow(raw), ncol=ncol(raw))

for(i in 1:nrow(raw)){
  for(j in 1:ncol(raw)){
    ind <- which(conds==conds[j])
    resids[i,j] <- exprs(raw)[i,j]-mean(exprs(raw)[i,ind])
  }
}

iND <- which(featureCategory(raw)=="Undetermined", arr.ind=TRUE)
iD <- which(featureCategory(raw)!="Undetermined", arr.ind=TRUE)
boxes_sub50 <- list("observed"=-resids[iD], "non-detects"=-resids[iND])

png(file="Boxplots_resids.png", 
    res = 300,
    width = 12, 
    height = 12, 
    units = "cm")

boxplot_sub50 <- boxplot(boxes_sub50, 
                         ylim=c(-5,5),
                         main="Comparison of residuals before imputation",
                         ylab=expression(paste("Ct residuals")))

dev.off()

#impute values using Single
imputed_single <- qpcrImpute(raw, 
                             dj = NULL, 
                             pyfit = NULL, 
                             groupVars="sampleType",
                             batch = NULL, 
                             tol = 1, 
                             iterMax = 100, 
                             outform=c("Single"), 
                             formula = NULL,
                             linkglm = "logit")

#create boxplots with single imputed values
conds_imputed_single <- paste(pData(imputed_single)$sampleType)

resids_imputed_single <- matrix(nrow=nrow(imputed_single), ncol=ncol(imputed_single))

for(i in 1:nrow(imputed_single)){
  for(j in 1:ncol(imputed_single)){
    ind <- which(conds==conds_imputed_single[j])
    resids[i,j] <- exprs(imputed_single)[i,j]-mean(exprs(imputed_single)[i,ind])
  }
}

iI_imputed_single <- which(featureCategory(imputed_single)=="Imputed", arr.ind=TRUE)
iD_mputed_single <- which(featureCategory(imputed_single)!="Imputed", arr.ind=TRUE)
boxes_imputed_single <- list("observed"=-resids[iD_mputed_single], "imputed"=-resids[iI_imputed_single])


png(file="Boxplots_resids_imputed.png", 
    res = 300,
    width = 12, 
    height = 12, 
    units = "cm")

boxplot_imputed <- boxplot(boxes_imputed_single, 
                           ylim=c(-5,5),
                           main="Comparison of residuals after imputation",
                           ylab=expression(paste("Ct residuals")))
dev.off()
