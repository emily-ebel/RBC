######################################################
#######  NORMALIZE PARASITE DATA  ########
######################################################

library(data.table) 
library(plyr)

# Read in data
# PMR = parasite multiplication rate      b=growth, a=invasion
# 'norm' is a simple normalization to 1111 from a given week, used later; 'raw' is the straight PMR
# Infected_0th = Day 0 parasitemia;   Infected_1st = Day 1 parasitemia
# raw1111x = raw 1111 PMR that week 
setwd("/Users/Emily/Desktop/Files")
pf1 <- read.table("3D7_data_1.14.20.csv", header = TRUE, sep = ",")  # I think this is the data with W9/10/11 normalization
pf1 <- pf1[,c(1,3,5,7,9,11,12,13,14,15)]
names(pf1) <- c("Sample","PMR_a_raw.3D7","PMR_b_raw.3D7","PMR_a_norm.3D7","PMR_b_norm.3D7","Infected_0th.3D7","Infected_1st.3D7","Week.3D7","raw1111a.3D7","raw1111b.3D7")       
pf2 <- read.table("TH_data_1.14.20.csv", header = TRUE, sep = ",") 
pf2 <- pf2[,c(1,3,5,7,9,11,12,13,14,15)]
names(pf2) <- c("Sample","PMR_a_raw.TH","PMR_b_raw.TH","PMR_a_norm.TH","PMR_b_norm.TH","Infected_0th.TH","Infected_1st.TH","Week.TH","raw1111a.TH","raw1111b.TH")       
data<-merge(pf1,pf2,all=TRUE)

# remove samples with no exome data/bone marrow transplant/no batch-effect data 
data <- subset(data, data$Sample != "6449KD" & data$Sample != "8715" & data$Sample != 'HETG')
data$Sample <-as.character(data$Sample)

#####################################################
############ LINEAR  BATCH EFFECT MODELS #############
#####################################################
# regress out effects of significant batch variables
# extract the remaining residuals -- these will be arithmetically corrected to finally represent fitness

## Add a variable for timing of experiment in weeks: 0= ran week collected, 1= ran next week, 2= ran 2 weeks later
# this is because parasite experiments failed in W9 (for 3D7) or W9 and W10 
data$RBC_age.3D7 = 0 
data$RBC_age.TH = 0 
oldw9 = c('3738','5691','6342','6449GR','6922')
oldw10 = c('3645','5609')
data <- within(data, RBC_age.TH[Sample %in% oldw9] <- 2)
data <- within(data, RBC_age.TH[Sample %in% oldw10] <- 2)
data <- within(data, RBC_age.3D7[Sample %in% oldw9] <- 1)

# ## Also add a variable for samples from W1-3, where PMRb comes from D3 -> D5 instead of D1 -> D3
early <- subset(data,data$Week.3D7<=2)
early <- subset(early,early$Week.3D7>0)
data$explen <- 3
data <- within(data, explen[Sample %in% early$Sample] <- 5)

# experimenter variable is NS for growth 


############
# growth
############
#### 3D7
#mod <- lm(data$PMR_b_raw.3D7~data$raw1111b.3D7+data$RBC_age.3D7+data$Infected_1st.3D7+data$explen+data$Week.3D7, na.action = na.exclude)
# removed NS variables to get
mod <- lm(data$PMR_b_raw.3D7~data$raw1111b.3D7+data$RBC_age.3D7, na.action = na.exclude)
summary(mod) #R2 0.3984 
data$res3D7b <- residuals(mod) # get residuals -- these are the measures of parasite growth after the above batch effects are removed
#summary(lm(data$res3D7b~data$Week.3D7)) # Week NS after other variables regressed out 

#### TH.026.09
#mod <- lm(data$PMR_b_raw.TH~data$raw1111b.TH+data$RBC_age.TH+data$Infected_1st.TH+data$explen+data$Week.TH, na.action = na.exclude)
mod <- lm(data$PMR_b_raw.TH~data$raw1111b.TH, na.action = na.exclude)
summary(mod)
data$resTHb <- residuals(mod) # get residuals 
#summary(lm(data$resTHb~data$Week.TH)) # Week NS after other variables regressed out 


############
# invasion
############
datainv <- subset(data,is.na(data$PMR_a_raw.3D7)==FALSE) # not all samples have invasion data
# try adding an experimenter variable
datainv$experimenter <- '0' # Emily
temp <- subset(datainv,datainv$Week.3D7>=26); temp2 <- subset(datainv,datainv$Week.3D7 <26)
temp$experimenter <- '1' #Bikash
datainv<-rbind(temp,temp2)

#### 3D7
mod <- lm(datainv$PMR_a_raw.3D7~datainv$raw1111a.3D7+datainv$Infected_0th.3D7+datainv$experimenter, na.action = na.exclude)
summary(mod)
datainv$res3D7a <- residuals(mod) # get residuals 

#### TH
mod <- lm(datainv$PMR_a_raw.TH~datainv$raw1111a.TH+datainv$Infected_0th.TH+datainv$RBC_age.TH, na.action = na.exclude)
summary(mod)
datainv$resTHa <- residuals(mod) # get residuals 



#####################################################
#############  ARITHMETIC CORRECTION - PT 1 ##############
# Residuals can be negative or positive (with mean 0), but we'd prefer to have them as percentages 
# To do this, we'll create a linear model to convert residuals to percentages based on a few set points 
#####################################################

### GROWTH
# 1111 is fixed at 100% (for now)
# The 5 HbAS donors + the HE donor are the most extreme: we'll set their percentages relative to 1111 from the respective week (aka, PMR_b_norm.XX from above)

datatrain <- subset(data,data$Sample %in% c("3185","3711","4612","7596","7995",'5678','1111_W4'))
# 3D7
datax3D7b <- datatrain$res3D7b
datay3D7b <- datatrain$PMR_b_norm.3D7
percent3D7bmod <- lm(datay3D7b~datax3D7b) 
new.df <- data.frame(datax3D7b=data$res3D7b)
data$percent3D7b <- predict(percent3D7bmod,new.df)*100
# TH
dataxTHb <- datatrain$resTHb
datayTHb <- datatrain$PMR_b_norm.TH
percentTHbmod <- lm(datayTHb~dataxTHb) 
new.df <- data.frame(dataxTHb=data$resTHb)
data$percentTHb <- predict(percentTHbmod,new.df)*100

### INVASION
# 1111 is fixed at 100% (for now)
# The G6PD_severe + the HE donor are the most extreme: we'll set their percentages relative to 1111 from the respective week (aka, PMR_b_norm.XX from above)
datainvtrain <- subset(datainv,datainv$Sample %in% c("6479",'5678','1111_W4'))
# 3D7
datainvx3D7a <- datainvtrain$res3D7a
datainvy3D7a <- datainvtrain$PMR_a_norm.3D7
percent3D7amod <- lm(datainvy3D7a~datainvx3D7a) 
new.df <- data.frame(datainvx3D7a=datainv$res3D7a)
datainv$percent3D7a <- predict(percent3D7amod,new.df)*100
# TH
datainvxTHa <- datainvtrain$resTHa
datainvyTHa <- datainvtrain$PMR_a_norm.TH
percentTHamod <- lm(datainvyTHa~datainvxTHa) 
new.df <- data.frame(datainvxTHa=datainv$resTHa)
datainv$percentTHa <- predict(percentTHamod,new.df)*100



#####################################################
#############  ARITHMETIC CORRECTION - PT 2 ##############
# We want to set the median non-carrier value to 100%
# Before doing this, remove repeated samples 
# and divide samples into carrier type 
#####################################################

#### Remove repeated samples, keeping the 1 with the most phenotypic data  
data_norep <- subset(data,data$Sample != "1111_W1" & data$Sample != "1111_W2" & data$Sample != "1111_W3" & data$Sample != "1111_W5" & data$Sample != "1111_W6" & data$Sample != "1111_W7" & data$Sample != "1111_W8" & data$Sample != "1111_W9" & data$Sample != "1111_W10" & data$Sample !=  "1111_W11" & data$Sample != "1111_W12" & data$Sample !=  "1111_W13" & data$Sample !=  "1111_W14" & data$Sample !=  "1111_W15" & data$Sample !=  "1111_W16" & data$Sample !=  "1111_W17" & data$Sample !=  "1111_W18" & data$Sample !=  "1111_W19" & data$Sample !=  "1111_W20" & data$Sample !=  "1111_W21" & data$Sample !=  "1111_W22" & data$Sample !=  "1111_W23" & data$Sample !=  "1111_W24" & data$Sample !=  "1111_W25" & data$Sample !=  "1111_W26" & data$Sample !=  "1111_W27")
data_norep <- subset(data_norep,data_norep$Sample != "2222_W17" & data_norep$Sample != "2222_W18" & data_norep$Sample != "2222_W25" & data_norep$Sample != "3333_W17" & data_norep$Sample != "4278_W2")
data_norep <- subset(data_norep,data_norep$Sample != "3730_W26" & data_norep$Sample != "3804_W27" & data_norep$Sample != "4278_W26" & data_norep$Sample != "5083_W26" & data_norep$Sample != "6443_W26" & data_norep$Sample != "7160_W27" & data_norep$Sample != "7496_W26" & data_norep$Sample != "8597_W26" & data_norep$Sample != '9172_W26')

#### Split samples by carrier type 
# Add the carrier info from data_norm.csv
sample_info <-  read.table("data_norm.csv", header = TRUE, sep = ","); sample_info=sample_info[,1:2] 
data_norep <- merge(data_norep,sample_info,by="Sample")

#### Divide samples by carrier status 
datahbas <- subset(data_norep,data_norep$Carrier =='HbAS')
dataHE <- subset(data_norep,data_norep$Carrier =='HE')
datanc <- subset(data_norep,data_norep$Carrier == 'no')
datahbac <- subset(data_norep,data_norep$Carrier =='HbAC')
dataalphahet <- subset(data_norep,data_norep$Carrier =='alphahet')
dataalphahom <- subset(data_norep,data_norep$Carrier =='alphahom')
datag6pdlow <- subset(data_norep,data_norep$Carrier == 'g6pdmild')
datag6pdhigh <- subset(data_norep,data_norep$Carrier =='g6pdsevere')

#### Standardize GROWTH so that healthy mean is 100%
datahbas$percent3D7b <- datahbas$percent3D7b/mean(datanc$percent3D7b)*100
dataHE$percent3D7b <- dataHE$percent3D7b/mean(datanc$percent3D7b)*100
datahbac$percent3D7b <- datahbac$percent3D7b/mean(datanc$percent3D7b)*100
dataalphahet$percent3D7b <- dataalphahet$percent3D7b/mean(datanc$percent3D7b)*100
dataalphahom$percent3D7b <- dataalphahom$percent3D7b/mean(datanc$percent3D7b)*100
datag6pdlow$percent3D7b <- datag6pdlow$percent3D7b/mean(datanc$percent3D7b)*100
datag6pdhigh$percent3D7b <- datag6pdhigh$percent3D7b/mean(datanc$percent3D7b)*100
datanc$percent3D7b <- datanc$percent3D7b/mean(datanc$percent3D7b)*100
#
datahbas$percentTHb <- datahbas$percentTHb/mean(datanc$percentTHb)*100
dataHE$percentTHb <- dataHE$percentTHb/mean(datanc$percentTHb)*100
datahbac$percentTHb <- datahbac$percentTHb/mean(datanc$percentTHb)*100
dataalphahet$percentTHb <- dataalphahet$percentTHb/mean(datanc$percentTHb)*100
dataalphahom$percentTHb <- dataalphahom$percentTHb/mean(datanc$percentTHb)*100
datag6pdlow$percentTHb <- datag6pdlow$percentTHb/mean(datanc$percentTHb)*100
datag6pdhigh$percentTHb <- datag6pdhigh$percentTHb/mean(datanc$percentTHb)*100
datanc$percentTHb <- datanc$percentTHb/mean(datanc$percentTHb)*100


## Now do it for invasion
# remove repeats 
datainv_norep <- subset(datainv,datainv$Sample != "1111_W1" & datainv$Sample != "1111_W2" & datainv$Sample != "1111_W3" & datainv$Sample != "1111_W5" & datainv$Sample != "1111_W6" & datainv$Sample != "1111_W7" & datainv$Sample != "1111_W8" & datainv$Sample != "1111_W9" & datainv$Sample != "1111_W10" & datainv$Sample !=  "1111_W11" & datainv$Sample != "1111_W12" & datainv$Sample !=  "1111_W13" & datainv$Sample !=  "1111_W14" & datainv$Sample !=  "1111_W15" & datainv$Sample !=  "1111_W16" & datainv$Sample !=  "1111_W17" & datainv$Sample !=  "1111_W18" & datainv$Sample !=  "1111_W19" & datainv$Sample !=  "1111_W20" & datainv$Sample !=  "1111_W21" & datainv$Sample !=  "1111_W22" & datainv$Sample !=  "1111_W23" & datainv$Sample !=  "1111_W24" & datainv$Sample !=  "1111_W25" & datainv$Sample !=  "1111_W26" & datainv$Sample !=  "1111_W27")
datainv_norep <- subset(datainv_norep,datainv_norep$Sample != "2222_W17" & datainv_norep$Sample != "2222_W18" & datainv_norep$Sample != "2222_W25" & datainv_norep$Sample != "3333_W17" & datainv_norep$Sample != "4278_W2")
datainv_norep <- subset(datainv_norep,datainv_norep$Sample != "3730_W26" & datainv_norep$Sample != "3804_W27" & datainv_norep$Sample != "4278_W26" & datainv_norep$Sample != "5083_W26" & datainv_norep$Sample != "6443_W26" & datainv_norep$Sample != "7160_W27" & datainv_norep$Sample != "7496_W26" & datainv_norep$Sample != "8597_W26" & datainv_norep$Sample != '9172_W26')
#### Divide samples by carrier status 
datainv_norep <- merge(datainv_norep,sample_info,by="Sample")
datainvhbas <- subset(datainv_norep,datainv_norep$Carrier =='HbAS')
datainvHE <- subset(datainv_norep,datainv_norep$Carrier =='HE')
datainvnc <- subset(datainv_norep,datainv_norep$Carrier == 'no')
datainvhbac <- subset(datainv_norep,datainv_norep$Carrier =='HbAC')
datainvalphahet <- subset(datainv_norep,datainv_norep$Carrier =='alphahet')
datainvalphahom <- subset(datainv_norep,datainv_norep$Carrier =='alphahom')
datainvg6pdlow <- subset(datainv_norep,datainv_norep$Carrier == 'g6pdmild')
datainvg6pdhigh <- subset(datainv_norep,datainv_norep$Carrier =='g6pdsevere')
#
datainvhbas$percent3D7a <- datainvhbas$percent3D7a/mean(datainvnc$percent3D7a)*100
datainvHE$percent3D7a <- datainvHE$percent3D7a/mean(datainvnc$percent3D7a)*100
datainvhbac$percent3D7a <- datainvhbac$percent3D7a/mean(datainvnc$percent3D7a)*100
datainvalphahet$percent3D7a <- datainvalphahet$percent3D7a/mean(datainvnc$percent3D7a)*100
datainvalphahom$percent3D7a <- datainvalphahom$percent3D7a/mean(datainvnc$percent3D7a)*100
datainvg6pdlow$percent3D7a <- datainvg6pdlow$percent3D7a/mean(datainvnc$percent3D7a)*100
datainvg6pdhigh$percent3D7a <- datainvg6pdhigh$percent3D7a/mean(datainvnc$percent3D7a)*100
datainvnc$percent3D7a <- datainvnc$percent3D7a/mean(datainvnc$percent3D7a)*100
#
datainvhbas$percentTHa <- datainvhbas$percentTHa/mean(datainvnc$percentTHa)*100
datainvHE$percentTHa <- datainvHE$percentTHa/mean(datainvnc$percentTHa)*100
datainvhbac$percentTHa <- datainvhbac$percentTHa/mean(datainvnc$percentTHa)*100
datainvalphahet$percentTHa <- datainvalphahet$percentTHa/mean(datainvnc$percentTHa)*100
datainvalphahom$percentTHa <- datainvalphahom$percentTHa/mean(datainvnc$percentTHa)*100
datainvg6pdlow$percentTHa <- datainvg6pdlow$percentTHa/mean(datainvnc$percentTHa)*100
datainvg6pdhigh$percentTHa <- datainvg6pdhigh$percentTHa/mean(datainvnc$percentTHa)*100
datainvnc$percentTHa <- datainvnc$percentTHa/mean(datainvnc$percentTHa)*100

#### SAVE (FOR LASSO)
growth_df <- as.data.frame(cbind(datanc$Sample,datanc$percent3D7b,datanc$percentTHb)); names(growth_df ) <- c("Sample","percent3D7b","percentTHb")
invasion_df <- as.data.frame(cbind(datainvnc$Sample,datainvnc$percent3D7a,datainvnc$percentTHa)); names(invasion_df) <- c("Sample","percent3D7a","percentTHa")
# fix incorrect sample names to be compatible with other scripts
growth_df$Sample <- gsub('1987DA', '1987W5', growth_df$Sample); growth_df$Sample <- gsub('1987ME', '1987W7', growth_df$Sample); growth_df$Sample <- gsub('6420', '5420', growth_df$Sample)
invasion_df$Sample <- gsub('1987DA', '1987W5', invasion_df$Sample); invasion_df$Sample <- gsub('1987ME', '1987W7', invasion_df$Sample); invasion_df$Sample <- gsub('6420', '5420', invasion_df$Sample)

write.table(growth_df,file="growth_df.csv",sep=",",col.names=TRUE,row.names=FALSE)
write.table(invasion_df,file="invasion_df.csv",sep=",",col.names=TRUE,row.names=FALSE)

