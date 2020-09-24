
##Load data 
data<-read.csv("2020_09_23_CO2_QAQC.csv")

##Format data
#Constants
R= 0.08205601 # gas constant in L*atm/mol*K
dens=0.9998395 # density of freshwater

#Input variables
temp_equil=data$temp_equil        # temperature of water immediately after equilibriaum in C
temp_samp=data$temp_samp         # temperature of water at the time of sampling in C
press_samp=data$press_samp        # pressure at the time of sampling in atm
vol_hs=data$vol_hs            # volume of headspace in mL
vol_samp=data$vol_samp          # volume of water sample in mL
CO2_pre=data$CO2_pre           # pCO2 of headspace before equilibrium (zero if using zero air) # uatm or ppmv 
CO2_post=data$CO2_post          # pCO2 of hs after equilibrium  # uatm or ppmv
A=data$alkalinity                 # alkalinity in umol/L

Calculate_CO2<-function(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post, A){
  StmpCO2 <-StmCO2fromSamp(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post)
  StmpCO2.umol <-Fw(tempC, StmpCO2)
  K1<-K1calc(temp_equil)
  K2<-K2calc(temp_equil)
  Carb1<-Carbfrom_C_A(K1, K2, StmpCO2.umol , A)
  delta_DIC<-DIC_correction(CO2_pre,StmpCO2,vol_hs,temp_equil,vol_samp)
  DIC_corr<-Carb1$D+delta_DIC
  Carb2<-Carbfrom_D_A(K1, K2, DIC_corr, A)
}


Calculate_CO2(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post,A)

############
## Export
############
write.csv(sum_file_final, "2020_09_23_CO2.csv",row.names = FALSE)



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
  StmpCO2 <- (CO2_post*KH.equil+(hsRatio*(CO2_post-CO2_pre)/molV))/KH.samp
  }

##### FUNCTION TO CONVERT pCO2 from uatm to umol/L #####  
Fw <- function(tempC, StmpCO2){
  temp_equil.K <- temp_equil + 273.15
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
Carbfrom_C_A <- function(K1, K2, StmpCO2.umol , A){
  H <- (((-K1*StmpCO2.umol ))-sqrt(((K1*StmpCO2.umol )^2)-(4*-1*A*2*K1*K2*StmpCO2.umol )))/(2*-1*A)
  pH <- -1*log10(H)
  B <- (K1*StmpCO2.umol )/H
  Ca <- (K2*B)/H
  D <- StmpCO2.umol  + B + Ca
  A_check <- B + 2*Ca
  Carb1 <- list(H, pH, StmpCO2.umol , B, Ca, D, A, A_check)
  names(Carb1) <- c("H", "pH", "StmpCO2.umol ", "B", "Ca", "D", "A","A check")
  Carb1
  }

### Function to calculate correction of perturbed sample for DIC before equilibration
## Originally from Dickenson et al. 2007--Chapter 4, highlighted again by Koschorreck et al. for freshwater

DIC_correction<-function(CO2_pre,CO2_post,vol_hs,temp_equil,vol_samp){
  delta_co2<-(((CO2_pre/1000/1000)-(CO2_post/1000/1000))*(vol_hs/1000))/R*(temp_equil + 273.15)
  delta_DIC<-delta_co2/(vol_samp*dens)
  delta_DIC
  }

### Apply correction to DIC data from Carb1####

DIC_corr<-DIC$D+delta_DIC


Carbfrom_D_A <- function(K1, K2, DIC_corr, A){
  a <- A
  b <- K1*(A-DIC_corr)
  c <- (A-(2*DIC_corr))*K1*K2
  H_t <- ((-1*b)+sqrt((b^2)-(4*a*c)))/(2*a)
  
  pH_t <- -1*log10(H_t)
  
  B_t <- (DIC_corr*K1*H_t)/((H_t^2)+(K1*H_t)+(K1*K2))
  Ca_t <- (DIC_corr*K1*K2)/((H_t^2)+(K1*H_t)+(K1*K2))
  
  C_t <- (H_t*B_t)/K1
  D2 <- C_t + B_t + Ca_t
  A2_check <- B_t + 2*Ca_t
  D_check <- C_t + B_t + Ca_t
  Carb2 <- list(H_t, pH_t, B_t, Ca_t, C_t, D_check, A2_check)
  names(Carb2) <- c("H", "pH", "B", "Ca", "CO2", "D_check", "A check")
  return(Carb2)
  }

#################################################
####Knit functions together to calculate CO2######
#################################################


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





