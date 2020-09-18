## Nyack Data Summary
library(plyr)
library(dplyr)
library(readr)     
library(ggplot2)
library(cowplot)
library(reshape2)
library(lubridate)
library(xts)
library(dygraphs)
library(tidyr)
library(purrr)
library(imputeTS)

######################
## Import and clean
######################
## setwd to where downloaded files are - eventually set to each individual location on SensorOne network
## Note: All files have different header and different columns so need to import individually

MET <- read.delim("CR1000_HA07_Met.dat", skip=1,sep=",")
MET <- MET[-c(1:2),]

#############################
## Extract DateTime & Light
############################
MET$DateTime <- as.POSIXct(as.character(MET$TIMESTAMP), format="%Y-%m-%d %H:%M:%S")
MET <- MET[,c("DateTime","PAR_Den")]
MET$DateTime <- force_tz(MET$DateTime, 'US/Mountain')

MET <- subset(MET, DateTime >= as.POSIXct("2018-08-15 00:00:00") & DateTime <= as.POSIXct("2018-09-06 00:00:00"))
## Weird middle of the night signals
MET$PAR_Den <- as.numeric(as.character(MET$PAR_Den))
MET$PAR_Den <- ifelse(MET$PAR_Den < 1, yes=0, no=MET$PAR_Den)

######################
## Spread time data
######################
## Create a dataset with the appropriate time values that you want
spread_ts = function(x, samp_freq){
  re = regexec("([0-9]+)([A-Z])",samp_freq)[[1]]
  if(-1%in%re){
    stop("Please enter a correct string")
  }else{
    ml = attr(re,"match.length")
    nn = as.numeric(substr(samp_freq, re[2], ml[2]))
    uu = substr(samp_freq, re[3], ml[1])
    if(uu=="D"){td = 24*60*60*nn
    }else if(uu=="H"){td = 60*60*nn
    }else if(uu=="M"){td = 60*nn
    }else if(uu=="S"){td = nn
    }else{stop("Please enter a correct string")}
  }
  seq(x[1], x[length(x)], by=td)
}

## Want 5 minute data
sampfreq = "5M"
# Generate a full timeseries at a given interval
DateTimeNew <- spread_ts(MET$DateTime, sampfreq)
xspread <- tibble(DateTime=DateTimeNew)
# Create a new table with NAs for missing time steps
Rnew <- merge(MET,xspread,by="DateTime", all.y=TRUE)
## Use Approx
Rnew$PAR_filled <- na.approx(Rnew$PAR)

## Visualize
ggplot(Rnew, aes(DateTime, PAR_filled))+geom_point()


#############################
## Combine with DO data
###########################
DO <- read.csv("2018_09_04_Beaver_DOcompiled.csv", header=T)
DO$DateTime <- as.POSIXct(as.character(DO$Time), format="%Y-%m-%d %H:%M:%S")

## Add one minute to light times
Rnew$DateTime <- Rnew$DateTime + (60)

## Merge
DOL <- merge(DO,Rnew,by="DateTime")
names(DOL)
DOL <- DOL[,c("Time","Temp","DO","PAR_filled")]

ggplot(DOL, aes(Time, PAR_filled))+geom_point()

write.csv(DOL, "2018_09_04_Beaver_DOLight.csv")

