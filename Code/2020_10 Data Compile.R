## Compiling discrete and continuous data from 2018 diel sampling
## JRB

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "readxl"), require, character.only=T)

## Reset wd
setwd("../Data")

################
## Import Data
################

## Import QAQC'd headspace equilibration data
dat <- ldply(list.files(pattern = "*_QAQC.csv"), function(filename) {
  d <- read.csv(filename)
  d$file <- filename
  return(d)
})
names(dat)
# subset to relevant columns
hs <- dat[,c("SampleID","CO2","WaterTemp_C","Baro_inHg")]
# import alkalinity and pH data
A_pH <- read.csv("Diel_pH_Alk_mean_reps.csv", header=T) ## mean of triplicate values for Oak
# convert from mg L-1 CaCO3 to mol m-3 CaCO3
A_pH$Alk_mol_m3 <- A_pH$Alk_mgLCaCO3/(60.008*1000)
# merge with iso data
hs <- left_join(hs, A_pH, by="SampleID")
## Separate SampleID to create DateTime
sapply(hs, class)
hs$SampleID <- as.character(hs$SampleID)
hs <- hs %>% 
  separate(SampleID, c("Site","Date","Time"))

hs$Date <- revalue(hs$Date, c("822" = "2018-08-22","823"="2018-08-23","830" = "2018-08-30","831"="2018-08-31",
                                "926" = "2018-09-26","927"="2018-09-27","928"="2018-09-28"))
hs$Time <- as.numeric(hs$Time)
hs$Time <- ifelse(hs$Time < 1000, yes= paste(substring(hs$Time, 1, 1), substring(hs$Time, 2, 3), sep = ":"),
                   no = ifelse(hs$Time >= 1000, yes=paste(substring(hs$Time, 1, 2), substring(hs$Time, 3, 4), sep = ":"), no=NA))
hs$Time <- revalue(hs$Time, c("0:"="0:00"))

hs$DateTime <- lubridate::ymd_hm(paste(hs$Date,hs$Time))
hs$DateTime <- as.POSIXct(as.character(hs$DateTime), format="%Y-%m-%d %H:%M:%S",tz = "UTC")
hs <- hs[order(hs$Site, hs$Date, hs$Time),]


## Import DO data
DOdat <- ldply(list.files(pattern = "*_DOcompiled.csv"), function(filename) {
  d <- read.csv(filename)
  d$file <- filename
  return(d)
})
names(DOdat)
# subset data
DOdat <- DOdat[,c("Time","Temp","DO","file")]
colnames(DOdat) <- c("DateTime","temp","oxy","file")
# split and process
sDO <- split(DOdat, DOdat$file)
names(sDO) <- c("Blaine","Beaver","OakLolomai","OakWillow")
# recombine
sDO <- lapply(sDO, function(x) return(x[,c("DateTime","temp","oxy")]))
dDO <- ldply(sDO, data.frame)
dDO$DateTime <- as.POSIXct(as.character(dDO$DateTime), format = "%Y-%m-%d %H:%M:%S")
## Round DateTime to the nearest 5 minute interval to help with integration
dDO$DateTime <- floor_date(dDO$DateTime, "5 minutes")
# rename columns
colnames(dDO) <- c("Site","DateTime","WaterTempC","DO_mgL")


##########
## Merge
##########
m <- merge(hs, dDO, by=c("Site","DateTime"))
## some DO/water temp data missing despite presence in dDO? something is wrong


