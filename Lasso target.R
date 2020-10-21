######################################################
#####  PERFORM LASSO ON NON-CARRIER PHENOTYPES  ####
#####  and/or NON-CARRIER GENOTYPES in 23 GENES  ####
######################################################

library(glmnet)
library(dplyr)
library(MLmetrics)
library(data.table)
library(MASS)

###########################################################
#####################  STEP 1   ###########################
###########  PREPARE DATA MATRIX FOR LASSO   ##############
###########################################################

#### Read in phenotype data:
setwd("/Users/Emily/Desktop/Files")
pheno0 <- read.table("data_norm.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE) 
# add osmo50
Osmo50 <-read.table("osmo50_norm.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE) 
pheno <- merge(pheno0,Osmo50,all=TRUE)
# get non-carriers only
pheno <- subset(pheno,pheno$Carrier=='no')
# Remove repeats 
pheno <- subset(pheno,pheno$Sample != "1111_W1" & pheno$Sample != "1111_W2" & pheno$Sample != "1111_W3" & pheno$Sample != "1111_W5" & pheno$Sample != "1111_W6" & pheno$Sample != "1111_W7" & pheno$Sample != "1111_W8" & pheno$Sample != "1111_W9" & pheno$Sample != "1111_W10" & pheno$Sample !=  "1111_W11" & pheno$Sample != "1111_W12" & pheno$Sample !=  "1111_W13" & pheno$Sample !=  "1111_W14" & pheno$Sample !=  "1111_W15" & pheno$Sample !=  "1111_W16" & pheno$Sample !=  "1111_W17" & pheno$Sample !=  "1111_W18" & pheno$Sample !=  "1111_W19" & pheno$Sample !=  "1111_W20" & pheno$Sample !=  "1111_W21" & pheno$Sample !=  "1111_W22" & pheno$Sample !=  "1111_W23" & pheno$Sample !=  "1111_W24" & pheno$Sample !=  "1111_W25" & pheno$Sample != "2222_W17" & pheno$Sample != "2222_W18" & pheno$Sample != '2222_W25' & pheno$Sample != "3333_W17" & pheno$Sample != "4278_W2" & pheno$Sample != "1111_W26" & pheno$Sample != "1111_W27" & pheno$Sample != "3730_W26" & pheno$Sample != "3804_W27" & pheno$Sample != "7496_W26"& pheno$Sample != "5083_W26"& pheno$Sample != "9172_W26"& pheno$Sample != "6449KD")
# drop non-phenotype, non-variable, or incomplete columns, plus Matrix columns where associations are driven by 1 outlier
pheno = pheno[ , -which(names(pheno) %in% c("Week.Advia","Repeat","Carrier","Matrix_9","Matrix_1","Matrix_3","Matrix_2","Retic_percent","Retic_num"))]
# Standardize sample names, Matrix 5, and make Sample an integer 
pheno$Sample <- gsub('1111_W4', '1111', pheno$Sample); pheno$Sample <- gsub('2222_W20', '2222', pheno$Sample); pheno$Sample <- gsub('3333_W18', '3333', pheno$Sample); pheno$Sample <- gsub('4278_W18', '4278', pheno$Sample); pheno$Sample <- gsub('1987DA', '1987W5', pheno$Sample); pheno$Sample <- gsub('1987ME', '1987W7', pheno$Sample); pheno$Sample <- gsub('6420', '5420', pheno$Sample)
pheno$Matrix_5 = pheno$Matrix_5/max(pheno$Matrix_5)*100
pheno$Sample <- as.integer(pheno$Sample)

#### Read in curated genotype data with parasite data:
# But first, choose whether you want to model invasion (a) or growth (b). 
#The sample sizes are different, so the pruning is different

geno_inv <- read.table("invasion_geno.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
geno_growth <- read.table("growth_geno.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE) 
# These dfs contain the top 10 PCs and unlinked SNPs/indels from 23 genes for each sample
# SNPs/indels are given shorthand names of 'v' for variant + number. 
# Actual info (e.g. position) on all variants in target genes  is available in varinfo.csv. Variants from this file have been pruned for linkage (see Methods) to yield the files above 

#### Combine pheno+geno for either invasion or growth
data_table_1 = data.table(pheno, key="Sample")
#data_table_2 = data.table(geno_inv, key="Sample")
data_table_2 = data.table(geno_growth, key="Sample")
x2 <- merge(data_table_1, data_table_2)
x <- as.data.frame(x2) 
# x contains all the data to perform Lasso, but it can be subsetted to test different things 

#### Select phenotypes only 
x.test <- x[,2:26]
#### Select genotypes and phenotypes 
x.test <- x[,2:(dim(x)[[2]]-2)] #all

#x.test = x.df3 
#### Create a copy of the chosen variables (this is important because permutations will shuffle the data later)
x.overall <- x.test; overall$Sample <- x$Sample

#### PICK THE PARASITE FEATURE YOU WANT TO TEST ####
#x.overall$parasite <- x$percent3D7a
#x.overall$parasite <- x$percentTHa
x.overall$parasite <- x$percent3D7b
#x.overall$parasite <- x$percentTHb

#### Prepare matrix for model 
x.mat = data.matrix(x.test[,1:dim( x.test)[[2]]]) 
# this is actual input into Lasso







###########################################################
#####################  STEP 2   ###########################
###############  RUN LASSO 1000 TIMES   ###################
###########################################################

# Create empty vectors to save MSE. R2, and variables with non-zero coefficients from each iteration
real.mse.vec = c()
real.r2.vec = c() # real.r2.1.vec = c() 
coef.nonzero.list = list()

# Now run Lasso 1000 times to predict parasite fitness (x.overall$parasite) from your chosen variable in x.mat
for (i in 1:1000){
  cvfit = cv.glmnet(x.mat, x.overall$parasite, standardize=TRUE) ## ACTUAL MODEL, fit for 100 lambdas (automatic cross-validation in choice of lambda)
  fitted = predict(cvfit, newx = x.mat, s='lambda.min') # use this model to make predictions of parasite fitness 
  plot(fitted,x.overall$parasite,xlab="Predicted parasite",ylab="Actual parasite",main=paste(i)) # look at how good the predictions are 
  # find and save the MSE for this iteration 
  real.mse = mean((fitted - x.overall$parasite)^2) 
  real.mse.vec = c(real.mse.vec,real.mse)
  # find and save the r2 for this iteration 
  testmod<-lm(x.overall$parasite~fitted)
  real.r2 = summary(testmod)$adj.r.squared 
  if (is.na(coef(testmod)[[2]])==TRUE){real.r2=0} else if (coef(testmod)[[2]]<0){real.r2=-real.r2}
  real.r2.vec = c(real.r2.vec,real.r2)
  # pull out variables with non-zero coefficients 
  coef = as.data.frame(as.matrix(coef(cvfit, s = "lambda.min"))); names(coef) = c('beta')
  coef.nonzero = subset(coef,coef$beta!=0) 
  coef.nonzero.list[[i]] = cbind(rownames(coef.nonzero),coef.nonzero$beta)  #### WE NEED THIS 
}


#### Summarize how often each x variable was selected in this set of 1000 models
num.coef.real.vec = c() 
names.coef.real.vec = c()
for (elem in coef.nonzero.list){
  num.coef.real.vec = c(num.coef.real.vec,dim(elem)[[1]]) 
  names.coef.real.vec = c(names.coef.real.vec,elem[,1]) 
}
tab = as.data.frame(table(names.coef.real.vec)); names(tab) <-c("Predictor","Selected")
# tab: contains the name of each predictor + the number of times it was chosen 

#### Find the median R2/MSE out of the 1000 models 
real.mse = median(real.mse.vec)
real.r2 = median(real.r2.vec)







###########################################################
#####################  STEP 3   ###########################
####  PERMUTE PARASITE DATA TO ASSESS SIGNIFICANCE   ######
###########################################################

num_perm = 1000 
num_cv = 10  # cross-validations per perm

# vectors to store summary stats for each permutation 
mse.perm = c()
r2.perm= c()

for (p in 1:num_perm){
  # permute parasite data 
  perm.parasite = sample(x.overall$parasite, replace=FALSE)
  # vectors to store summary stats for each CV within this permutation 
  mse.perm.cv = c()
  r2.perm.cv = c()
  # do the CV
  for (i in 1:num_cv) {
    cvfit = cv.glmnet(x.mat, perm.parasite, standardize=TRUE) 
    fitted = predict(cvfit, newx = x.mat, s='lambda.min') 
    # save summary stats from this cv 
    mse = mean((fitted - perm.parasite)^2)
    mse.perm.cv = c(mse.perm.cv, mse)
    testmod<-lm(perm.parasite~fitted); r2.perm.cv = c(r2.perm.cv, summary(testmod)$adj.r.squared)
  }
  # if any r2 from the CV of this perm are NA (meaning no model), set them to 0
  for (elem in 1:length(r2.perm.cv)){if (is.na(r2.perm.cv[[elem]])==TRUE){r2.perm.cv[[elem]]=0}}
  # save summary stats (median of CVs of this perm)
  mse.perm = c(mse.perm, median(mse.perm.cv))
  r2.perm= c(r2.perm, median(r2.perm.cv))
}


###################

##### Compare real and permuted medians 
# MSE
msepval=sum(mse.perm<=real.mse,na.rm=TRUE)/length(mse.perm)
msepval
hist(mse.perm,xlim=c(0,max(mse.perm)),main=paste("median MSE",msepval),xlab="")
abline(v=real.mse,lty=2)

## R2
r2pval = sum(r2.perm>=real.r2,na.rm=TRUE)/length(r2.perm)
r2pval
hist(r2.perm,cex.axis=1.4,col='gray',cex.lab=1.5,xlab = expression(paste('R'^"2")), ylab="", main="",bty = 'n',xlim=c(0,1))
abline(v=real.r2,lty=2,col=rgb(245, 179, 66, max = 255, alpha = 190, names="orange"),lwd=3)
mtext(text="Permutations",side=2,line=2.5,cex=1.5)#side 1 = bottom