
library(plyr)
library(dplyr)
library(readr)     
library(ggplot2)
library(cowplot)
library(reshape2)
library(xts)
library(dygraphs)
library(tidyr)

###########################################
##Load Data
###########################################

## Prior to running the rest of the code
## 1) Set the wd to folder location of the raw Picarro dat files and meta data
## 2) Specify meta file name
meta_file_name <- "2020_08_04_13CO2.test2.txt"
## 3) Specify start time of Picarro run (using Picarro time)
picarro_start_time <- "2020-08-04 19:00:00"
## 4) Specify notes file name
notes_file_name <- "Chamber_test.run_08.04.2020.csv"

###########################################
##Calculate pCO2 of water from pCO2 headspace
###########################################

##### SOLUBILITY CONSTANT FOR CO2 (KH) #####
# Use Henry's Law and parameters from Weiss 1974 to calculate KH in units of mol L-1 atm-1 as a function of temperature (also in Demarty et al 2011; many others)
KH.CO2 <- function(temp_equil){
  tempK <- temp_equil + 273.15
  KH.CO2 <- exp( -58.0931 + (90.5069*(100/tempK)) +22.2940*log((tempK/100), base=exp(1)) )
  KH.CO2
}

##### FUNCTION TO ESTIMATE STREAM pCO2 from headspace CO2 #####
# R= gas constant (0.08205601) in L*atm/mol*K
# temp_equil= temperature of water immediately after equilibriaum in C
# temp_samp= temperature of water at the time of sampling in C
# press_samp= pressure at the time of sampling in atm
# vol_hs= volume of headspace in mL
# vol_samp= volume of water sample in mL
# CO2_pre=pCO2 of headspace before equilibrium (zero if using zero air) # uatm or ppmv 
# CO2_post=pCO2 of hs after equilibrium  # uatm or ppmv
StmCO2fromSamp <- function(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post){
  temp_equil.K <- temp_equil + 273.15
  molV <- 0.08205601*(temp_equil.K)*(press_samp) # L mol-1
  hsRatio <- vol_hs/vol_samp
  KH.equil <- KH.CO2(temp_equil) # mol L-1 atm-1
  KH.samp <- KH.CO2(temp_samp) # mol L-1 atm-1
  StmCO2 <- (CO2_post*KH.equil+(hsRatio*(CO2_post-CO2_pre)/molV))/KH.samp
  StmpCO2
}

##### FUNCTION TO CONVERT pCO2 from uatm to umol/L #####  
Fw <- function(tempC, Cw.uatm){
  StmpCO2.umol <- StmpCO2*(1/0.08205601)*(1/(10^-3))*(1/temp_equil.K)/1000
  StmpCO2.umol
}





############
## Export
###########
write.csv(sum_file_final, "2020_08_04_13CO2_QAQC.csv",row.names = FALSE)



