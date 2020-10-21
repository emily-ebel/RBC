######################################################
###  OSMOTIC FRAGILITY SUMMARY & NORMALIZATION  ######
######################################################
# Clean up and fit a sigmoid model to the raw data
# Then normalize model curves by variation in 1111 (who has OF data every week)

setwd("/Users/Emily/Desktop/Files")
raw <-read.table("raw_OF_data.csv", header = TRUE, sep = ",") 
# raw: each row is a sample; columns indicate osmolality (2 reps per); cells have absorption (OD540) data
# outlier points were manually converted to NA previously


#### Convert OD540 to hemoglobin content (g/dL, from lysis)
hb <- raw #make a copy
hb[3:34] <- lapply(hb[3:34], function(x) x*200/100/68)  ### conversion 
hb <- hb[-c(3,19,18,34) ] # Remove tubes 1 (highest [NaCl]) and 16 (lowest [NaCl]) since these were dropped in most weeks
# normalize hemoglobin content by max per row/sample
maxes<-apply(hb[3:30],1,max, na.rm=TRUE) #max of each row
hb_norm0 <- sweep(hb[3:30], 1, maxes, '/') # apply to data cols
# hb_norm0 contains the normalized hgb data only 
hb_norm <- cbind(hb[1:2],hb_norm0)  
# hb_norm has 2 header cols added


#### Fit sigmoid models to the hemoglobin data: basically we want to build a finer-scale model with more x-values 

# this vector has the x-values we measured: osmolality values for each tube (columns in hb_norm0)
x <- c(245.5,210,196,182,168,154,147,140,133,126,112,105,91,84,245.5,210,196,182,168,154,147,140,133,126,112,105,91,84)
# sample plot of hemolysis by osmolality
y<- hb_norm[hb_norm$Sample == '1111_W12', ] 
plot(x,y[3:30],xlab="Osmolality",ylab="% Hemolysis")
# here are the x-values we want to model
x_mod <- seq(80, 246, length = 500) # we will estimate a y for each of these x
# for both x vectors, normalize osmolality to "relative tonicity" (where 100 is the max theoretical salt solution, 308 mosmol)
x <- x/308*100
x_mod <- x_mod/308*100 #max  called for 

# Function to determine sigmoid 
sigmoid <- function(col) {
  df1 = data.frame(v1=x,v2=col)
  names(df1) = c("Osmolality","Hemolysis")
  summary(mod1<-nls(Hemolysis ~ SSlogis( log(Osmolality), Asym, xmid, scal),data=df1))
  z<-predict(mod1, data.frame(Osmolality = x_mod))
  # scale the predicted y-values (z) so max is 100
  z = z/max(z)
  # peek at the data if you want
  #plot(x_mod,z)
 # return(z)
}

# Apply this function to data from all samples 
sig0 <-apply(hb_norm0,1,sigmoid) #Pass each row (sample) of data to sigmoid function; in output q, Each column is a sample (with names from hb[1:1]), with 500 rows corresponding to x_mod
# clean up 
sig <- t(sig0)
sig <- sig*100
sig_nometa <- sig
sig <- cbind(hb[1:2],sig)
sig <- as.data.frame(sig) 
## sig has the UNNORMALIZED model curves. Each row is a sample; each column corresponds to an x_mod (relative tonicity) value; cells are hemolysis

# Take a look at first 10 curves if you want
# plot(x_mod,sig[2,3:502],col="white",xlab="Relative tonicity",ylab="hemolysis")
# for (i in 1:10){
#   lines(x_mod,sig[i,3:502])
# }


############################################
###### Find the 'osmo50' of each sample ######
# i.e., the relative tonicity (x-value) corresponding to 50% hemolysis

# empty matrix to store osmo50 data 
osmo50 <- as.data.frame(matrix(, nrow = 1, ncol = dim(sig_nometa)[[1]]))

# calculate osmo50 for each sample
for (i in 1:dim(sig_nometa)[[1]]) { 
  x=sig_nometa[i,] # the whole vector of hemolysis across 500 x-values 
  o50y= 50 # half of max hemolysis for the sample 
  o50z=which(abs(x-o50y)==min(abs(x-o50y))) #index of the x(hemolysis) value closest to the osmo50 (half of max hemolysis)
  osmo50[1,i] <- x_mod[o50z]
}

# clean up
names(osmo50) <- sig$Sample
osmo50t<-t(osmo50)
osmo50_unnorm <- cbind(rownames(osmo50t), data.frame(osmo50t, row.names=NULL))
colnames(osmo50_unnorm) <- c("Sample",'Osmo50') 
# osmo50_final has UNNORMALIZED osmo50 for each sample 

# save these unnormalized values 
setwd("/Users/Emily/Desktop/Files")
write.table(osmo50_unnorm,file="osmo50_unnorm.csv",sep=",",col.names=TRUE,row.names=FALSE)



#### #### #### #### #### #### #### 
#### NORMALIZATION STRATEGY #####
#### #### #### #### #### #### #### 

# 1. Find 1111 average osmo50
# 2. Calculate 1111's deviation each week from the average 
# 3. Apply -(deviation) to the curves for 1111 for each week
# 4. Plot to make sure it works
# 5. Apply to all samples
# 6. Plot again
# 7. Recalculate Osmo50


# 1. Find 1111 average osmo50
library(data.table)
osmo50_unnorm_1111 <- osmo50_unnorm[osmo50_unnorm$Sample %like% "1111", ]
osmo50_1111 <- mean(osmo50_unnorm_1111$Osmo50, na.rm=TRUE)  
addw3 <- data.frame('1111_W3',osmo50_1111) # this week something was wrong with 1111 OF data, so just make it the average (need something for normalization)
names(addw3)<-names(osmo50_unnorm_1111 )
osmo50_unnorm_1111  <- rbind(osmo50_unnorm_1111 ,addw3)

# 2. Calculate the deviation each week from the average 
osmo50_unnorm_1111$diff <- osmo50_unnorm_1111$Osmo50-osmo50_1111
osmo50_unnorm_1111$corrected <- osmo50_unnorm_1111$Osmo50-osmo50_unnorm_1111$diff
osmo50_unnorm_1111$Week <- c (1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,3)

# 3. Apply -(deviation) to the curves for 1111 for each week
# We need to adjust the x-values (x_mod), not the y-values (hemolysis)
x_mod_df_1111 <- matrix(, nrow = 24, ncol = 500) # each row holds the x values for sig's y values (which will be modified by normalization)
for(row in 1:24){
  x_mod_df_1111[row, ] <- x_mod # original x-values for 1111
}
for (row in 1:dim(x_mod_df_1111)[[1]]){
  x_mod_df_1111[row,] <- x_mod_df_1111[row,]-osmo50_unnorm_1111[row,3] # adjustd x-values for 1111
}

# 5. Apply the deviation to all samples
x_mod_df <- matrix(, nrow = 149, ncol = 500)
# pre-fill the matrix
for(row in 1:149){
  x_mod_df[row, ] <- x_mod
}  
x_mod_df<-cbind(sig[,1:2],x_mod_df)
# add the deviation 
for (i in 1:dim(x_mod_df)[[1]]){
  osmo50_unnorm_1111row = subset(osmo50_unnorm_1111,osmo50_unnorm_1111$Week==x_mod_df[i,2])
  x_mod_df[i,3:502] <- as.numeric(x_mod_df[i,3:502])-as.numeric(osmo50_unnorm_1111row[3])
}  



# 4. Plot to make sure it works (1111)
sig1111<- sig[sig$Sample %like% "1111", ]
# #before
plot(0, ylab="Hemolysis",main='',xlim=c(30,100), ylim=c(0,100), xlab='Relative Tonicity',cex.axis=1.6,cex.lab=1.8, bty = 'n', pch = '')
for (i in 1:dim(sig1111)[[1]]){
  lines(x_mod,sig1111[i,3:502])
}
abline(h=50)
abline(v=osmo50_1111)
#after
plot(0, ylab="Hemolysis",main='',xlim=c(30,100), ylim=c(0,100), xlab='Relative Tonicity',cex.axis=1.6,cex.lab=1.8, bty = 'n', pch = '')
for (i in 1:dim(sig1111)[[1]]){
  lines(x_mod_df_1111[i,],sig1111[i,3:502])}
abline(h=50)
abline(v=osmo50_1111)


# It works! So the normalized data curves are in x_mod_df (1 row per sample) and sig (1 row per sample)     


# 7. Recalculate osmo50 
# It works as long as x_mod[z] comes from x_mod_df[i,][z]
sig_nometa <- sig[3:502]
x_mod_df_nometa <- x_mod_df[3:502]
val50 <- sig_nometa[,1]/2 #target values -- half of max hemolysis (from row [new:column] 1)

#matrix to store data 
osmo50_norm <- as.data.frame(matrix(, nrow = 1, ncol = length(val50)))

for (i in 1:length(val50)) { # for every sample 
  x=sig_nometa[i,] # the whole vector of hemolysis 
  y=as.numeric(val50[i]) # half of max hemolysis for the sample 
  z=which(abs(x-y)==min(abs(x-y))) #index of the x(hemolysis) value closest to the osmo50 (half of max hemolysis)
  # osmo50[1,i] = sig_nometa[,z]$xval
  osmo50_norm[1,i] <- x_mod_df_nometa[i,z]
}

names(osmo50_norm) <- sig$Sample
m<-t(osmo50_norm)
osmo50_norm <- cbind(rownames(m), data.frame(m, row.names=NULL))
colnames(osmo50_norm) <- c("Sample",'Osmo50') #### 'nn' has NORMALIZED osmo50 for each sample 

# save normalized osmo50
write.table(osmo50_norm,file="osmo50_norm.csv",sep=",",col.names=TRUE,row.names=FALSE)
