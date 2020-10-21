######################################################
#######  PLOT EKTA DATA  ########
######################################################

library(TSdist)

# Raw data format: two columns per sample: one is osmolality (=mm NaCl*2) and the other is the deformability index
setwd("/Users/Emily/Desktop/Files")
d<-read.table("raw_ekta_data.csv", header = TRUE, sep = ",") 
d[d == 0] <- NA


######################################################
########## FUNCTION TO SMOOTH RAW DATA ###############
# good for a quick glance at the data
# but does not take into account week-to-week batch effects 
######################################################

myfunction <- function(sample){ # pass the sample label  
  # pull out all the data
  y_idx = match(paste("X",sample,sep=''),names(d)) # find the index of this sample (column with y (DI) values)
  x_idx = y_idx - 1 # index of the column with x (osmolality) value  
  data = as.data.frame(cbind( d[,x_idx] , d[,y_idx] ))
  data <- data[complete.cases(data), ]
  names(data) = c("x","y")
  # Some sample curves are odd on the edges and need slightly different treatment -- determined from looking at individual graphs of raw data
  need_lower = c('3121','2751','3711','7700','6922','9743','7677','6443') 
  hightail1 <- c('1287'); hightail2 <- c('2112')  

    if (sample %in% need_lower) { data1 <- subset(data,data$x>80 & data$x<490)
    } else if (sample %in% hightail1) { data1 <- subset(data,data$x>110 & data$x<470) 
  } else if (sample %in% hightail2) { data1 <- subset(data,data$x>110 & data$x<430) 
  } else {data1 <- subset(data,data$x>110 & data$x<490) 
  }

  # build the model
  model <- lm(data1$y ~ poly(data1$x,20)) # use 20 degress to get smoothest fit 
  predicted.intervals <- as.data.frame(predict(model,data.frame(x=data1$x),interval='confidence', level=0.99))
  # return the model
  model_df <- unique(as.data.frame(cbind(data1$x,predicted.intervals$fit)))
  return(model_df)
}

# example use
z <- myfunction("1111_W1")
plot(z$V1,z$V2,xlab="Osmolality (mm NaCl*2)", ylab="Deformability Index")

######################################################

######################################################
########## FUNCTION TO GET COORDS OF 3 KEY POINTS ####
# from the raw data 
# these are used to split the curve up into chunks that are normalized separately
#####################################################

myfunction2 <- function(curve,samp){ # pass a df with x-y points of the curve (output of previous function)
  # max points - based on max y value 
  DImax <- subset(curve,curve$y==max(curve$y))$y
  Omax <- subset(curve,curve$y==max(curve$y))$x
  
  # min points -- start by basing on min value in first half of curve 
  # BUT if the min value is the first value, do an inflection function
  exceptions = c('6922','9743','7677') # these guys need an inflection function 
  if ((samp %in% exceptions)==FALSE){
    data3 <- subset(curve,curve$x < 300)
    DImin <- subset(data3,data3$y==min(data3$y))$y
    Omin <- subset(data3,data3$y==min(data3$y))$x
  } else {
    # we need some other way of calculating Omin
    if (samp == '9743'){DImin=0.115;Omin=117} # just eyeballed these since there are only 3
    if (samp == '6922'){DImin=0.158;Omin=132}
    if (samp == '7677'){DImin=0.1225;Omin=115}
  }
  
  # hyper points
  if (samp != '3874') {
    data4 <- subset(curve,curve$x > Omax) 
    rowindex <- which(abs(data4$y-(DImax/2))==min(abs(data4$y-(DImax/2))))
    DIhyper <- data4[rowindex,]$y
    Ohyper <- data4[rowindex,]$x 
  } else {Ohyper=420;DIhyper=DImax/2} # 3874 has cutoff data 
  
  
  return(c(DImax,Omax,DImin,Omin,DIhyper,Ohyper))
}



################################################################
###############  PLOT NORMALIZED CURVES ##########################
################################################################
## --> shifts the curves to correct for week-to-week variation 


### Read in factors that correspond to vertical and horizontal correction for each of 3 points, based on normalization factors by week from median_normalization.R
factors <- read.table("ektanormfactors.csv", header = TRUE, sep = ",") 
names(factors)[1] = 'Week.Advia'
## Also read in which week each sample was collected in
sw <- read.table("sampleweek.csv", header = TRUE, sep = ",") 

### Read in a data frame that has each donor's carrier status 
data <- read.table("data_norm.csv", header = TRUE, sep = ",")  
data$Sample <- as.character(data$Sample)
# remove repeats, keeping the one with the most phenotypic data (1111_W4, 2222_W20, 3333_W18, 4278_W18)
data2 <- subset(data,data$Sample != "1111_W1" & data$Sample != "1111_W2" & data$Sample != "1111_W3" & data$Sample != "1111_W5" & data$Sample != "1111_W6" & data$Sample != "1111_W7" & data$Sample != "1111_W8" & data$Sample != "1111_W9" & data$Sample != "1111_W10" & data$Sample !=  "1111_W11" & data$Sample != "1111_W12" & data$Sample !=  "1111_W13" & data$Sample !=  "1111_W14" & data$Sample !=  "1111_W15" & data$Sample !=  "1111_W16" & data$Sample !=  "1111_W17" & data$Sample !=  "1111_W18" & data$Sample !=  "1111_W19" & data$Sample !=  "1111_W20" & data$Sample !=  "1111_W21" & data$Sample !=  "1111_W22" & data$Sample !=  "1111_W23" & data$Sample !=  "1111_W24" & data$Sample !=  "1111_W25" & data$Sample !=  "1111_W26" & data$Sample !=  "1111_W27")
data2 <- subset(data2,data2$Sample != "2222_W17" & data2$Sample != "2222_W18" & data2$Sample != "2222_W25" & data2$Sample != "3333_W17" & data2$Sample != "4278_W2")
data2 <- subset(data2,data2$Sample != "3730_W26" & data2$Sample != "3804_W27" & data2$Sample != "4278_W26" & data2$Sample != "5083_W26" & data2$Sample != "6443_W26" & data2$Sample != "7160_W27" & data2$Sample != "7496_W26" & data2$Sample != "8597_W26" & data2$Sample != '9172_W26')
data2 <- subset(data2, data2$Sample != "6449KD") # remove the sample with missing exome data
# get vector of donors in each category
noncarrier <- subset(data2,data2$Carrier == 'no')$Sample
hbas <- subset(data2,data2$Carrier == 'HbAS')$Sample
hbac <- subset(data2,data2$Carrier == 'HbAC')$Sample
alphahet <- subset(data2,data2$Carrier == 'alphahet')$Sample
alphahom <- subset(data2,data2$Carrier == 'alphahom')$Sample
HE <- subset(data2,data2$Carrier == 'HE')$Sample
g6pdlow <- subset(data2,data2$Carrier == 'g6pdmild')$Sample
g6pdhigh <- subset(data2,data2$Carrier == 'g6pdsevere')$Sample



# set up empty plots 
plot(0,bty="n",xlab="",ylab="",cex=0.5,col="white",ylim=c(0,0.65),xlim=c(50,250),pch=21,cex.axis=1.2,cex.lab=1)
mtext(text="Deformability Index",side=2,line=2.5,cex=1.6)#side 1 = bottom
mtext(text="NaCl (mM)",side=1,line=2.5,cex=1.6)#side 1 = bottom

#### Go through the ekta data and plot the desired samples 
for (i in 1:((dim(d)[2])/2)) {
  # GET SAMPLE NAME
  col2 = i*2
  samp = names(d)[col2]
  if (substring(samp,1,1)=='X'){samp = substring(samp,2)}
  ### IS IT SOMETHING WE WANT TO PLOT? edit this as needed -- just one example shown
  if ( (samp %in% noncarrier) | (samp %in% hbac) ) {  
      # PASS NAME TO CURVE FUNCTION
      z <- myfunction(samp) # z just contains the points to plot
      names(z) <- c('x','y') 
      # PASS CURVE AND NAME TO POINTS FUNCTION
      points <- myfunction2(z,samp)
      
      #####################################
      # NOW NORMALIZE ALL POINTS IN CURVE #
      #####################################
    
      # make a copy of z
      newz <- z
      ## Get the factors for this sample
      week <- subset(sw,sw$Sample==samp); week = as.character(week$Week.Advia)
      normfacs <-subset(factors, factors$Week.Advia == week); normfacs <- normfacs[2:7]
      ## Consider each point -- this takes a couple seconds per sample
      for (row in 1:dim(z)[[1]]){
        xval = z[row,1]
        yval = z[row,2]
        # if point has x-value at or below Omin, apply Omin/DI min normalization
        if (xval <= points[4]){
              # apply Omin/DI min normalization
              newx <- xval*normfacs[4]
              newy <- yval*normfacs[3]
              newz[row,1] <- newx; newz[row,2] <- newy
          # else if point has x-value at or below Ohyper, apply Ohyper/DIhyper normalization
        } else if (xval >= points[6]){
              newx <- xval*normfacs[6]
              newy <- yval*normfacs[5]
              newz[row,1] <- newx; newz[row,2] <- newy
              # else if point has x-value between Omin and Omax (inclusive),
        } else if ((xval >= points[4]) & ( xval <= points[2]) ) {
              # calculate distance between each point
              distOmin <- EuclideanDistance(c(xval,yval), c(points[4],points[3]))
              distOmax <- EuclideanDistance(c(xval,yval), c(points[2],points[1]))
              fracOmin <- 1-(distOmin/(distOmin+distOmax)); fracOmax = 1-fracOmin
              # Weight the factors by distance
              xfactor = fracOmin*normfacs[4] + fracOmax*normfacs[2]
              yfactor = fracOmin*normfacs[3] + fracOmax*normfacs[1]
              # Apply that normalization
              newx <- xval*xfactor
              newy <- yval*yfactor
              newz[row,1] <- newx; newz[row,2] <- newy
          # if point has x-value between Omax and Ohyper (inclusive),
        } else if ((xval >= points[2]) &  (xval <= points[6])) {
              # calculate distance between each point
              distOmax <- EuclideanDistance(c(xval,yval), c(points[2],points[1]))
              distOhyper <- EuclideanDistance(c(xval,yval), c(points[6],points[5]))
              fracOmax <- 1-(distOmax/(distOmax+distOhyper)); fracOhyper = 1-fracOmax
              # Weight the factors by distance
              xfactor = fracOhyper*normfacs[6] + fracOmax*normfacs[2]
              yfactor = fracOhyper*normfacs[5] + fracOmax*normfacs[1]
              # Apply that normalization
              newx <- xval*xfactor
              newy <- yval*yfactor
              newz[row,1] <- newx; newz[row,2] <- newy }
      } # end 'consider each point' loop
    
      ## PLOT
      ## Convert x from Osmolality units to mM units (just divide by 2) 
      newz$x <- newz$x/2
       if (samp %in% noncarrier) { lines(newz$x,newz$y,col=rgb(142, 211, 245, max = 255),lwd=2) } # blue
      else if (samp %in% hbas){ lines(newz$x,newz$y,col=rgb(1,0,0,0.7), lwd=2) }
      else if (samp %in% HE){ lines(newz$x,newz$y,col=rgb(0, 0, 0, max = 255, alpha = 90, names = "gray"),lwd=2) }
      else if (samp %in% hbac){ lines(newz$x,newz$y,col=rgb(255, 177, 61, max = 255, alpha = 180, names = "orange"),lwd=2) }
      else if (samp %in% g6pdhigh){ lines(newz$x,newz$y,col=rgb(20, 97, 8, max = 255, alpha = 255, names = "green"),lwd=2,lty=3) }
      else if (samp %in% g6pdlow){  lines(newz$x,newz$y,col=rgb(38, 156, 20, max = 255, alpha = 180, names = "green"),lwd=2)  }
      else if (samp %in% alphahom){ lines(newz$x,newz$y,col=rgb(174, 97, 237, max = 255, alpha = 255, names = "purple"),lwd=2,lty=3) }
      else if (samp %in% alphahet){ lines(newz$x,newz$y,col=rgb(174, 97, 237, max = 255, alpha = 255, names = "purple"),lwd=2) }
}}
  
