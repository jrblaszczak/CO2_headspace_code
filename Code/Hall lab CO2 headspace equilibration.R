
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
Fw <- function(tempC, StmpCO2){
  StmpCO2.umol <- StmpCO2*(1/0.08205601)*(1/(10^-3))*(1/temp_equil.K)/1000
  StmpCO2.umol
}

###########################################
##Calculate DIC of water from pCO2 of water
###########################################

## Functions to determine equilibrium constants dependent on temperature of water when sampled
K1calc<- function(temp_samp) { 10^( (-3404.71/(273.15+temp_samp)) + 14.844 -0.033*(temp_samp+273.15) )}
K2calc<- function(temp_samp) { 10^( (-2902.39/(273.15+temp_samp)) + 6.498 -0.0238*(temp_samp+273.15) )}

### Function to solve carbonate chemistry using CO2 & Alkalinity
## See R Markdown file (XXXXXXX) for derivation of equations
# H=H+ 
# A=alkalinity in umol/L
# StmpCO2.umol=CO2 in umol/L 
# B=bicarbonate 
# Ca=carbonate 
# D=DIC in umol/L
# A_check=alkalinity check
Carbfrom_C_A <- function(K1, K2, StmpCO2.umol, A){
  H <- (((-K1*StmpCO2.umol))-sqrt(((K1*StmpCO2.umol)^2)-(4*-1*A*2*K1*K2*StmpCO2.umol)))/(2*-1*A)
  pH <- -1*log10(H)
  B <- (K1*StmpCO2.umol)/H
  Ca <- (K2*B)/H
  D <- StmpCO2.umol + B + Ca
  A_check <- B + 2*Ca
  
  l <- list(H_minus, pH, StmpCO2.umol, B, Ca, D, A, A_check)
  names(l) <- c("H", "pH", "StmpCO2.umol", "B", "Ca", "D", "A","A check")
  DIC<-D
  return(l)
}

### Function to calculate correction of perturbed sample for DIC before equilibration
## Originally from Dickenson et al. 2007--Chapter 4, highlighted again by Koschorreck et al. for freshwater

DIC_correction<-function(CO2_pre,CO2_post,vol_hs,temp_equil,vol_samp){
  delta_co2<-((CO2_pre-CO2_post)*vol_hs)/R*temp_equil
  delta_DIC<-delta_co2/(vol_samp*dens)
  delta_DIC
  }


Carbfrom_D_A <- function(K1, K2, DIC_corr, A){
  H <- ((-1*K1*(A-DIC_corr))+sqrt((B^2)-(4*A*(A-(2*DIC_corr))*K1*K2)))/(2*A)
  pH <- -1*log10(H)
  B <- (DIC_corr*K1*H)/((H^2)+(K1*H)+(K1*K2))
  Ca <- (DIC_corr*K1*K2)/((H^2)+(K1*H)+(K1*K2))
  CO2 <- (H*B)/K1
  D_check <- C + B + Ca
  A_check <- B + 2*Ca
  
  l <- list(H, pH, B, Ca, CO2, D_check, A_check)
  names(l) <- c("H", "pH", "B", "Ca", "C", "D_check", "A check")
  return(l)
  CO2_corr<-CO2
}

#################################################
####Knit functions together to calculate CO2######
#################################################
##Load data
#Constants
R= 0.08205601 # gas constant in L*atm/mol*K
dens=0.9998395 # density of freshwater

#Input variables
temp_equil=        # temperature of water immediately after equilibriaum in C
temp_samp=         # temperature of water at the time of sampling in C
press_samp=        # pressure at the time of sampling in atm
vol_hs=            # volume of headspace in mL
vol_samp=          # volume of water sample in mL
CO2_pre=           # pCO2 of headspace before equilibrium (zero if using zero air) # uatm or ppmv 
CO2_post=          # pCO2 of hs after equilibrium  # uatm or ppmv
A=                 # alkalinity in umol/L

Calculate_CO2<-function(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post, A){
  StmpCO2 <-StmCO2fromSamp(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post)
  StmpCO2.umol <-Fw(tempC, StmpCO2)
  K1<-K1calc(temp_equil)
  K2<-K2calc(temp_equil)
  DIC<-Carbfrom_C_A(K1, K2, StmpCO2.umol, A)
  delta_DIC<-DIC_correction(CO2_pre,CO2_post,vol_hs,temp_equil,vol_samp)
  DIC_corr<-DIC+delta_DIC
  CO2_corr<-Carbfrom_D_A(K1, K2, DIC_corr, A)
  }

Calculate_CO2(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post,A)





############
## Export
############
write.csv(sum_file_final, "2020_08_04_13CO2_QAQC.csv",row.names = FALSE)




