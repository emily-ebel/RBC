######################################################
#####  PERFORM LASSO ON NON-CARRIER PHENOTYPES  ####
#####  and/or NON-CARRIER GENOTYPES in 23 GENES  ####
######################################################

library(glmnet)
library(dplyr)
library(MLmetrics)
library(data.table)
library(MASS)
require(caret)

###########################################################
#####################  STEP 1   ###########################
###########  PREPARE DATA MATRIX FOR LASSO   ##############
###########################################################

#### Read in phenotype data:
setwd("/Users/Emily/Desktop/RBC-main")
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
# Standardize sample names, Matrix values, and make Sample an integer 
pheno$Matrix_sum = pheno$Matrix_4+pheno$Matrix_5+pheno$Matrix_6+pheno$Matrix_7+pheno$Matrix_8
pheno$Matrix_4 = pheno$Matrix_4/pheno$Matrix_sum*100
pheno$Matrix_5 = pheno$Matrix_5/pheno$Matrix_sum*100
pheno$Matrix_6 = pheno$Matrix_6/pheno$Matrix_sum*100
pheno$Matrix_7 = pheno$Matrix_7/pheno$Matrix_sum*100
pheno$Matrix_8 = pheno$Matrix_8/pheno$Matrix_sum*100
pheno<-pheno[,1:26] # drop the Matrix_sum column 
pheno$Sample <- gsub('1111_W4', '1111', pheno$Sample); pheno$Sample <- gsub('2222_W20', '2222', pheno$Sample); pheno$Sample <- gsub('3333_W18', '3333', pheno$Sample); pheno$Sample <- gsub('4278_W18', '4278', pheno$Sample); pheno$Sample <- gsub('1987DA', '1987W5', pheno$Sample); pheno$Sample <- gsub('1987ME', '1987W7', pheno$Sample); pheno$Sample <- gsub('6420', '5420', pheno$Sample)
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
#data_table_2 = data.table(geno_inv, key="Sample") ## uncomment this to test invasion 
data_table_2 = data.table(geno_growth, key="Sample") ## uncomment this to test growth
x2 <- merge(data_table_1, data_table_2)
x <- as.data.frame(x2) 
# remove the 5 siblings (leaving in the mother 2752)
fam = c('2814','3627','4379','6129','9387')
x<-subset(x,(x$Sample %in% fam)==FALSE)

# x now contains all the data to perform Lasso. 
# In the loops below, x will be split into train and test sets many times to perform k-folds CV.
# Also in the loops below, select which parasite phenotype to model!




###########################################################
#####################  STEP 2   ###########################
######  SPLIT DATA INTO TRAIN/TEST SETS MANY TIMES  #######
###########################################################


set.seed(123) 
kvalue = 10 # number of folds 
median.r2.train.kfold = c() # these vectors will store the overall results 
median.r2.test.kfold = c()
median.mse.train.kfold = c()
median.mse.test.kfold = c()
predictors_all = c()


for (split in 1:3) { #number of iterations, or splits of the data into folds 
  print(c("SPLIT",split))
  flds <- createFolds(x$Sample, k = kvalue, list = TRUE, returnTrain = FALSE)
  predictors_onesplit = c() # these vectors will store k results from this split of the data  
  median.r2.train = c()
  median.r2.test = c()
  median.mse.train = c()
  median.mse.test = c()
  
  #### For each fold,
  for (q in 1:kvalue){
    x.test = x[flds[[q]], ] # that fold is the test set 
    x.train = subset(x,(x$Sample %in% x.test$Sample == FALSE))  # and the other folds are the train set 
    
    #### CHOOSE A SET OF PREDICTORS TO TEST #####
    #### Select genotypes only
    #   x.select <- x.train[,27:(dim(x.train)[[2]]-2)] ## TRAIN
    #  x.select.test <- x.test[,27:(dim(x.test)[[2]]-2)] # TEST
    #### Select genotypes and phenotypes 
    x.select <- x.train[,2:(dim(x.train)[[2]]-2)] ## TRAIN
    x.select.test <- x.test[,2:(dim(x.test)[[2]]-2)] # TEST
    #### Select phenotypes only 
    # x.select <- x.train[,2:26] ## TRAIN
    #  x.select.test <- x.test[,2:26] # TEST
    
    ### PICK THE PARASITE PHENOTYPE: one of (percent3D7b, percentTHb) for growth or (percent3D7a, percentTHa) for invasion
    x.overall <- x.select; x.overall$Sample <- x.select$Sample; x.overall$parasite <- x.train$percent3D7b # TRAIN
    test.parasite <-x.test$percent3D7b # TEST 
    
    #### Prepare matrix for model - this is actual input into Lasso
    x.mat = data.matrix(x.select[,1:dim( x.select)[[2]]]) ## TRAIN 
    x.mat.test = data.matrix(x.select.test[,1:dim( x.select.test)[[2]]]) # TEST
    
    ###########################################################
    #####################  STEP 3   ###########################
    ##################  RUN LASSO    ##########################
    ###########################################################
    
    ### Lasso has an internal CV procedure to select lambda that relies on random splits of the input data (our train data)
    ### Since the splits are random, we'll run it 5 times and pick the median model out of 5

    real.mse.vec = c() # these vectors will store the 5 results from Lasso internal CV 
    real.r2.vec = c() 
    coef.nonzero.list = list()
    real.mse.vec.test = c()
    real.r2.vec.test = c() 
  
    for (i in 1:5){
      cvfit = cv.glmnet(x.mat, x.overall$parasite, standardize=TRUE) ## ACTUAL MODEL for TRAIN data, fit for 100 lambdas (automatic cross-validation in choice of lambda)
      fitted = predict(cvfit, newx = x.mat, s='lambda.min') # prediction accuracy for TRAIN data 
      fitted.test = predict(cvfit, newx = x.mat.test, s='lambda.min') # prediction accuracy for TEST data 
      # save MSE 
      real.mse = mean((fitted - x.overall$parasite)^2) # TRAIN
      real.mse.vec = c(real.mse.vec,real.mse)  # TRAIN
      real.mse.test = mean((fitted.test - test.parasite)^2) # TEST
      real.mse.vec.test = c(real.mse.vec.test,real.mse.test) # TEST
      # save r2 
      testmod<-lm(x.overall$parasite~fitted)  # TRAIN
      real.r2 = summary(testmod)$adj.r.squared   # TRAIN
      real.r2.vec = c(real.r2.vec,real.r2)  # TRAIN
      testmod.test<-lm(test.parasite~fitted.test) # TEST
      real.r2.test = summary(testmod.test)$adj.r.squared  # TEST
      real.r2.vec.test = c(real.r2.vec.test,real.r2.test) # TEST
      
      # pull out variables with non-zero coefficients (for TRAIN only)
      coef = as.data.frame(as.matrix(coef(cvfit, s = "lambda.min"))); names(coef) = c('beta')
      coef.nonzero = subset(coef,coef$beta!=0) 
      coef.nonzero.list[[i]] = cbind(rownames(coef.nonzero),coef.nonzero$beta)   
    }

    #### Save  the median R2/MSE from the 5 internal CV 
    real.mse = median(real.mse.vec) # TRAIN
    median.mse.train = c(median.mse.train,real.mse) # TRAIN
    real.r2 = median(real.r2.vec) # TRAIN
    median.r2.train = c(median.r2.train,real.r2) # TRAIN
    real.mse.test = median(real.mse.vec.test) # TEST
    real.r2.test = median(real.r2.vec.test) # TEST 
    median.r2.test = c(median.r2.test,real.r2.test) # TEST 
    median.mse.test = c(median.mse.test,real.mse.test) # TEST 
    
    ##### Save the predictors from the median model of the 5 internal CV  (TRAIN only)
    names.coef.real.vec = c()
    for (elem in coef.nonzero.list){
      names.coef.real.vec = c(names.coef.real.vec,elem[,1])
    }
    tab = as.data.frame(table(names.coef.real.vec)); names(tab) <-c("Predictor","Selected")
    predictors_onesplit_onesubset = as.character(tab$Predictor)
    predictors_onesplit = c(predictors_onesplit,predictors_onesplit_onesubset)
    
  } # end kvalue loop
  
  # save the mean MSE/R2 values from this set of 10 folds 
  median.r2.train.kfold = c(median.r2.train.kfold,mean(median.r2.train))
  median.r2.test.kfold = c(median.r2.test.kfold,mean(median.r2.test))
  median.mse.train.kfold = c(median.mse.train.kfold,mean(median.mse.train))
  median.mse.test.kfold = c(median.mse.test.kfold,mean(median.mse.test))
  # save the predictors 
  predictors_all=c(predictors_all,predictors_onesplit)
  
} #end split loop

# save the outputs so you don't have to run the code again 
df <- cbind(median.r2.train.kfold,median.r2.test.kfold,median.mse.train.kfold,median.mse.test.kfold)
df <- as.data.frame(df)
names(df)=c('r2.train','r2.test','mse.train','mse.test')
saveRDS(df, 'pheno.kmeans.stats.3D7b.rds') 
tab_cv = as.data.frame(table(predictors_all)); names(tab_cv) <-c("Predictor","Selected") # 'Selected' is the number of folds*number of splits in which the predictor was selected in the train model 
saveRDS(tab_cv, 'pheno.kmeans.predictors.3D7b.rds') 





###########################################################
#####################  STEP 4   ###########################
####  PERMUTE PARASITE DATA TO ASSESS SIGNIFICANCE   ######
###########################################################
## This code is the same as above, except for one line (222) and renamed variables 


require(caret)
set.seed(123)
num_perm = 1000
num_cv = 5  # internal cross-validations per permutation (5 above)
kvalue = 10
x.perm <- x # same data frame as above 

# vectors to store summary stats for each permutation 
perm.overall.r2.train.vec = c()
perm.overall.r2.test.vec = c()
perm.overall.mse.train.vec = c()
perm.overall.mse.test.vec = c()


for (p in 1:num_perm){
  print(c("PERM",p))
  ## THE LINE BELOW! PERMUTE PARASITE DATA by random sampling 
  x.perm$percentTHb = sample(x$percentTHb, replace=FALSE) ## ADJUST BASED ON PARASITE OUTCOME - ALSO ADJUST IN LASSO LOOP BELOW!                                                    
  median.r2.train.kfold.perm = c()                           
  median.r2.test.kfold.perm = c()
  median.mse.train.kfold.perm = c()
  median.mse.test.kfold.perm = c()
  
  for (split in 1:10) { #  
    print(c("SPLIT",split))
    flds <- createFolds(x.perm$Sample, k = kvalue, list = TRUE, returnTrain = FALSE)

    median.r2.train.perm = c()
    median.r2.test.perm = c()
    median.mse.train.perm = c()
    median.mse.test.perm = c()
    
    for (q in 1:kvalue){ # for each fold 
      x.test.perm = x.perm[flds[[q]], ]
      x.train.perm = subset(x.perm,(x.perm$Sample %in% x.test.perm$Sample == FALSE))  
      
      #### Select genotypes and phenotypes 
      #   x.select <- x.train[,2:(dim(x.train)[[2]]-2)] ## TRAIN
      #    x.select.test <- x.test[,2:(dim(x.test)[[2]]-2)] # TEST
      #### Select phenotypes only 
      x.select.train.perm <- x.train.perm[,2:26] ## TRAIN
      x.select.test.perm <- x.test.perm[,2:26] # TEST
      
      #### Prepare matrix for model - this is actual input into Lasso
      x.mat.train.perm = data.matrix(x.select.train.perm[,1:dim( x.select.train.perm)[[2]]]) 
      # but also do for test set 
      x.mat.test.perm = data.matrix(x.select.test.perm[,1:dim( x.select.test.perm)[[2]]]) 
      
      ### ACTUALLY RUN LASSO ####
      perm.mse.vec.train = c()
      perm.r2.vec.train = c() 
      # save summaries of test data (from train model)
      perm.mse.vec.test = c()
      perm.r2.vec.test = c() 
      
      for (i in 1:num_cv){
        cvfit = cv.glmnet(x.mat.train.perm, x.train.perm$percentTHb, standardize=TRUE) 
        fitted = predict(cvfit, newx = x.mat.train.perm, s='lambda.min')
        fitted.test = predict(cvfit, newx = x.mat.test.perm, s='lambda.min') 
        perm.mse.train = mean((fitted - x.train.perm$percentTHb)^2) 
        perm.mse.vec.train = c(perm.mse.vec.train,perm.mse.train)
        perm.mse.test = mean((fitted.test - x.test.perm$percentTHb)^2) 
        perm.mse.vec.test = c(perm.mse.vec.test,perm.mse.test)
        testmod<-lm(x.train.perm$percentTHb~fitted)
        perm.r2.train = summary(testmod)$adj.r.squared 
        perm.r2.vec.train = c(perm.r2.vec.train,perm.r2.train)
        testmod.test<-lm(x.test.perm$percentTHb~fitted.test)
        perm.r2.test = summary(testmod.test)$adj.r.squared 
        perm.r2.vec.test = c(perm.r2.vec.test,perm.r2.test)
      } # end internal Lasso CV loop
      
      perm.mse.train = median(perm.mse.vec.train)       # these vectors will be length k 
      perm.r2.train = median(perm.r2.vec.train)
      perm.mse.test = median(perm.mse.vec.test)
      perm.r2.test = median(perm.r2.vec.test)
      median.r2.train.perm = c(median.r2.train.perm,perm.r2.train)
      median.r2.test.perm = c(median.r2.test.perm,perm.r2.test)
      median.mse.train.perm = c(median.mse.train.perm,perm.mse.train)
      median.mse.test.perm = c(median.mse.test.perm,perm.mse.test)

    } # end kvalue loop
    
    # save the values from this split: these vectors will be length(splits) 
    median.r2.train.kfold.perm = c(median.r2.train.kfold.perm,mean(median.r2.train.perm))
    median.r2.test.kfold.perm = c(median.r2.test.kfold.perm,mean(median.r2.test.perm))
    median.mse.train.kfold.perm = c(median.mse.train.kfold.perm,mean(median.mse.train.perm))
    median.mse.test.kfold.perm = c(median.mse.test.kfold.perm,mean(median.mse.test.perm))     
  } # end split loop
  
  # Save an overall value from this permutation 
  perm.overall.r2.train = mean(median.r2.train.kfold.perm)
  perm.overall.r2.test = mean(median.r2.test.kfold.perm)
  perm.overall.mse.train = mean(median.mse.train.kfold.perm)
  perm.overall.mse.test = mean(median.mse.test.kfold.perm)
  perm.overall.r2.train.vec = c(perm.overall.r2.train.vec,perm.overall.r2.train)
  perm.overall.r2.test.vec = c(perm.overall.r2.test.vec,perm.overall.r2.test)
  perm.overall.mse.train.vec = c(perm.overall.mse.train.vec,perm.overall.mse.train)
  perm.overall.mse.test.vec = c(perm.overall.mse.test.vec,perm.overall.mse.test)
  
} # end permutation loop 

# save the vectors to access later
df <- cbind(perm.overall.r2.train.vec,perm.overall.r2.test.vec,perm.overall.mse.train.vec,perm.overall.mse.test.vec)
df <- as.data.frame(df)
names(df)=c('r2.train','r2.test','mse.train','mse.test')
saveRDS(df, 'pheno.kmeans.PERM.stats.3D7b.rds') 


