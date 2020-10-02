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
hs <- dat[,c("SampleID","CO2","WaterTemp_C","H2O_mL","air_mL","Baro_inHg")]
# import alkalinity and pH data
A_pH <- read.csv("Diel_pH_Alk_mean_reps.csv", header=T) ## mean of triplicate values for Oak
# merge with iso data
hs <- left_join(hs, A_pH, by="SampleID")
## Separate SampleID to create DateTime
sapply(hs, class)
#create a second sample ID to save
hs$Sample.ID <- hs$SampleID
#convert sample ID into site, date, time
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


#####################################################
## Format data to match requirements of Koschorreck
#####################################################

# If supplying a data frame, you can build it importing of a csv file
# Example: dataset <- read.csv("R_test_data.csv")
# The first row of this file must contain column names, then one row for each sample to be solved.
# The columns names must be:
# 
# 1. Sample.ID 
# 2. HS.pCO2.before
# 3. HS.pCO2.after
# 4. Temp.insitu
# 5. Temp.equil
# 6. Alkalinity.measured
# 7. Volume.gas
# 8. Volume.water
# 9. Bar.pressure
# 10. Constants
# 11. Salinity

names(hs)
# create each column if does not already exist
# hs$Sample.ID already exists
hs$HS.pCO2.before <- 0
hs$HS.pCO2.after <- hs$CO2 ## check if units correct
hs$Temp.insitu <- hs$WaterTemp_C
hs$Temp.equil <- hs$WaterTemp_C ## assuming the temperature is the same
hs$Alkalinity.measured <- hs$Alk_mgLCaCO3










