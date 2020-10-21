######################################################
#######  NORMALIZE BY MEDIAN OF NON-CARRIERS  ########
######################################################

# this applies to (1) all Advia data and (2) summary statistics from ekta 
# osmotic fragiliy and parasite data are normalized separately 
# plotting normalized ekta curves (in normalize_ekta.R) requires this script to be run first 

library(plyr)
library(dplyr)

# read in the unnormalized data
setwd("/Users/Emily/Desktop/Files")
data <-read.table("data_unnorm.csv", header = TRUE, sep = ",")

# separate out non-carriers
nc <- subset(data,data$Carrier == "no") 

# Collect the columns with data to be normalized
tonorm <- nc[,6:34]

# Find the median for each trait (across samples) in each week
traitmedians<-ddply(tonorm, .(Week.Advia), colwise(median,na.rm=TRUE))
traitmedians<-traitmedians[1:24,] # last row is NA b/c only parasite data were collected after W25
# traitmedians: each row is a week, each column a trait (except first column), each cell the median across samples for that week/trait

# Find the median value for each trait across weeks (this will be the new median PER week after normalization)
median_na <- function(x) median(x,na.rm=TRUE)
weeklymedians <- apply(traitmedians[,2:length(traitmedians)],2,median_na)
weeklymedians <- c(NA,weeklymedians) # first column is week -- treat median as NA 
weeklymedians <- as.data.frame(t(as.data.frame(weeklymedians)))
## weeklymedians: Each column is a trait, each cell is the median across weeks 

# drop any traits where the median is 0 (can't be normalized)
drop <- c() 
for (i in 2:length(weeklymedians)){
  if (weeklymedians[[i]]==0) { drop <- c(drop,names(weeklymedians)[[i]])}  }
traitmedians<-traitmedians[,-which(names(traitmedians) %in% drop)]
weeklymedians <-  as.vector(weeklymedians[,-which(names(weeklymedians) %in% drop)])

#### find a factor for each week/trait
# First create an empty df to save it in
Week.Advia <- as.vector(traitmedians$Week.Advia)
temp <- rep(0,length(Week.Advia)) # this is just a placeholder that will be deleted in a couple lines
facdf <- as.data.frame(cbind(Week.Advia,temp))
# find the factors
for (tr in 2:length(traitmedians)){ # for every trait 
  factors <- weeklymedians[[tr]]/traitmedians[,tr]  # each element is the factor to normalize for that week/trait
  facdf <- cbind(facdf,factors)
}  
# clean it up
facdf<-facdf[,-which(names(facdf) %in% c('temp'))]
names(facdf) <- names(traitmedians)
## facdf: rows are weeks, columns are traits, cells are the factor to multiply each sample by for that week/trait to normalize


#### Save the factors for the ekta points, which will be reused in plot_normalized_ekta.R
ektadf <- cbind(facdf$Week.Advia,facdf[,21:26])
write.table(ektadf,file="ektanormfactors.csv",sep=",",col.names=TRUE,row.names=FALSE)


#### now apply the factors to a copy of the data. Althougt the medians were found for non-carriers only, apply the normalization to everyone
# copy the data 
data2 <- data
# apply the normalization  
for (i in 1:(dim(data2)[[1]])){ # FOR EACH ROW OF DATA
    row = data2[i,]
    facrow <- subset(facdf,facdf$Week.Advia==row$Week.Advia) # Get the factors for the appropriate week 
    if (dim(facrow)[[1]] != 0) {
    for (tr in 2:length(facrow)){
      row.index <- match(names(facrow)[[tr]],names(row))
      if (facrow[[tr]]!=Inf){
        row[[row.index]] = row[[row.index]] * facrow[[tr]]}
      else { # if the factor is Inf, it means the value in traitmedians was 0 -- so instead of multiplying, just add the median 
        row[[row.index]] = row[[row.index]] + weeklymedians[[tr]] } }
    # SAVE
    data2[i,] = row } 
}


#### show that this makes the medians equal 
# before normalization (weekly median of non-carriers in red)
plot(nc$Week.Advia,nc$Ohyper)
temp <- nc %>% group_by(Week.Advia) %>% dplyr::summarize(Median = median(Ohyper, na.rm=TRUE))
points(temp$Week.Advia,temp$Median,col='red')
## after  (weekly median of non-carriers in blue)
plot(data2$Week.Advia,data2$Ohyper)
nc2 <- subset(data2,data2$Carrier == 'no')
temp <- nc2 %>% group_by(Week.Advia) %>% dplyr::summarize(Median = median(Ohyper, na.rm=TRUE))
points(temp$Week.Advia,temp$Median,col='blue')


#### save
write.table(data2,file="data_norm.csv",sep=",",col.names=TRUE,row.names=FALSE)
