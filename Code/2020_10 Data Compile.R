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
colnames(dDO) <- c()

