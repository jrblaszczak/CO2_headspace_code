
##Load data 
data<-read.csv("2020_10_02_HScorrection_preformat.csv")
data<-data[complete.cases(data),]
##Format data
#Constants
R= 0.08205601 # gas constant in L*atm/mol*K
dens=0.9998395 # density of freshwater


#Input variables
temp_equil=data$WaterTemp_C        # temperature of water immediately after equilibriaum in C
temp_samp=data$WaterTemp_C           # temperature of water at the time of sampling in C
press_samp=(data$Baro_inHg/29.921)        # pressure at the time of sampling in atm
vol_hs=data$air_mL            # volume of headspace in mL
vol_samp=data$H2O_mL         # volume of water sample in mL
CO2_pre=0          # pCO2 of headspace before equilibrium (zero if using zero air) # uatm or ppmv 
CO2_post=data$CO2          # pCO2 of hs after equilibrium  # uatm or ppmv
A=(data$Alk_mgLCaCO3/1000/100)                 # alkalinity in mol/L


Calculate_CO2<-function(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post, A){
 ##StmpCO2 <-StmCO2fromSamp(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post)
  StmpCO2.mol<-co2(temp_samp,CO2_post )
  K1<-K1calc(temp_equil)
  K2<-K2calc(temp_equil)
  Carb1<-Carbfrom_C_A(K1, K2, StmpCO2.mol , A)
  delta_DIC<-DIC_correction(CO2_pre,CO2_post,vol_hs,temp_equil,vol_samp)
  DIC_corr<-Carb1$D+delta_DIC
  Carb2<-Carbfrom_D_A(K1, K2, DIC_corr, A)
  Result <- list(Carb2$pH, StmpCO2.mol , delta_DIC,Carb2$CO2,Carb2$CO2/(KH.CO2(temp_equil)/1000000) )
  names(Result) <- c("pH", "StmpCO2_mol_original","delta_DIC","StmpCO2.mol_corrected","StmpCO2_corrected_ppmv")
  return(Result)
  }

##Returns a list with the corrected DIC (D) and CO2 
result<-Calculate_CO2(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post,A)
df<-as.data.frame(result)
df<-df[complete.cases(df),]

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
#
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
  molV <- R*(temp_equil.K)*(press_samp) # L mol-1
  hsRatio <- vol_hs/vol_samp
  KH.equil <- KH.CO2(temp_equil) # mol L-1 atm-1
  KH.samp <- KH.CO2(temp_samp) # mol L-1 atm-1
  StmpCO2 <- ((CO2_post*KH.equil)+(hsRatio*((CO2_post-CO2_pre)/molV)))/KH.samp
}

#####ppmv to mol/L
co2 <- function(temp_samp,CO2_post ){
      KH.samp <- KH.CO2(temp_samp)
      co2.mol.L<-  KH.samp*CO2_post/1000000
}
StmpCO2.mol<-co2(temp_samp,StmpCO2 )

###########################################
##Calculate DIC of water from pCO2 of water
###########################################

## Functions to determine equilibrium constants dependent on temperature of water when sampled
#From Milerno et al 2006 "Dissociation constants of carbonic acid in seawater as a function of salinity and temperature"
K1calc<- function(temp_equil) { 10^-(-126.34048+6320.813/(temp_equil+273.15)+19.568224*log(temp_equil+273.15))}
K2calc<- function(temp_equil) { 10^-(-90.18333+5143.692/(temp_equil+273.15)+14.613358*log(temp_equil+273.15))}


### Function to solve carbonate chemistry using CO2 & Alkalinity
## See R Markdown file (XXXXXXX) for derivation of equations
# H=H+ 
# A=alkalinity in mol/L
# StmpCO2.umol=CO2 in mol/L 
# B=bicarbonate 
# Ca=carbonate 
# D=DIC in umol/L
# A_check=alkalinity check
Carbfrom_C_A <- function(K1, K2, StmpCO2.mol , A){
  H <- (((-K1*StmpCO2.mol ))-sqrt(((K1*StmpCO2.mol )^2)-(4*-1*A*2*K1*K2*StmpCO2.mol )))/(2*-1*A)
  pH <- -1*log10((H))
  B <- (K1*StmpCO2.mol )/H
  Ca <- (K2*B)/H
  D <- StmpCO2.mol  + B + Ca
  A_check <- B + 2*Ca
  Carb1 <- list(H, pH, StmpCO2.mol , B, Ca, D, A, A_check)
  names(Carb1) <- c("H", "pH", "StmpCO2.mol ", "B", "Ca", "D", "A","A check")
  Carb1
  }

### Function to calculate correction of perturbed sample for DIC before equilibration
## Originally from Dickenson et al. 2007--Chapter 4, highlighted again by Koschorreck et al. for freshwater

DIC_correction<-function(CO2_pre,CO2_post,vol_hs,temp_equil,vol_samp){
  delta_DIC<-(((CO2_post-CO2_pre)/1000000))/(R*(temp_equil + 273.15))*(vol_hs/vol_samp)
  delta_DIC
  }

### Apply correction to DIC data from Carb1####

DIC_corr<-Carb1$D+delta_DIC


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
  DIC_corr<-Carb1$D+delta_DIC
  CO2_corr<-Carbfrom_D_A(K1, K2, DIC_corr, A)
  }


