## Test of Koschorrek HS equilibration correction

library(dplyr)
source("Koschorreck_Rheadspace.R")

## Import data and format again
setwd("../Data")
dat <- read.csv("2020_10_02_HScorrection_allformatted.csv")
dat <- dat[,-c(1)]
sapply(dat,class)
dat <- dat %>% mutate_if(is.integer,as.numeric)

## for now get rid of any NAs -- some alkalinity measurements missing because of titration error
dat <- na.omit(dat)

## Run code
pCO2 <- Rheadspace(dat)

## Export
write.csv(pCO2, "2020_10_02 Koschorrek hs equil correction results.csv")
