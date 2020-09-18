## Blaine Creek Diel

library(lubridate)
library(plyr)
library(dplyr)
library(readr)     
library(ggplot2)
library(cowplot)
library(reshape2)
library(data.table)
library(bayesbio)
library(tidyr)
library(tidyquant)

## Desired Outputs: 
## (1) a spreadsheet with DO, Temp, PAR, CO2, pH
## (2) a spreadsheet with discrete isotope timepoints

## Step 1: setwd

########################################################################################
## Import notes, CO2, O2, processed Picarro data, and pH sensor data
########################################################################################

###########
## Notes ##
###########

## merge again with other data if available
notes <- read.table("2018_08_23_BlaineDiel_notes.txt", sep=",",header=T)

notes <- notes %>% 
  separate(SampleID, c("Site","Date","Time"))

notes$Date <- revalue(notes$Date, c("822" = "2018-08-22","823"="2018-08-23"))
notes$Time <- as.numeric(notes$Time)
notes$Time <- ifelse(notes$Time < 1000, yes= paste(substring(notes$Time, 1, 1), substring(notes$Time, 2, 3), sep = ":"),
                     no = ifelse(notes$Time >= 1000, yes=paste(substring(notes$Time, 1, 2), substring(notes$Time, 3, 4), sep = ":"), no=NA))
notes$Time <- revalue(notes$Time, c("0:"="0:00"))

notes$DateTime <- lubridate::ymd_hm(paste(notes$Date,notes$Time))
notes$DateTime <- as.POSIXct(as.character(notes$DateTime), format="%Y-%m-%d %H:%M:%S")

########
## O2 ##
########
O2raw <- read.csv("2018_08_23_Blaine_DOcompiled.csv", header=T)
O2raw$DateTime <- as.POSIXct(as.character(O2raw$Time), format="%Y-%m-%d %H:%M:%S")

## Visualize
plot_grid(
  ggplot(O2raw, aes(DateTime, DO)) + geom_point(),
  ggplot(O2raw, aes(DateTime, Temp)) + geom_point(),
  ggplot(O2raw, aes(DateTime, Q)) + geom_point(),
  ncol=1, align="hv")

## Select Date Range of interest
O2 <- subset(O2raw, DateTime > "2018-08-22 05:30:00" & DateTime < "2018-08-23 06:30:00")

## Visualize
ggplot(O2, aes(DateTime, DO)) +
  geom_point()+
  geom_point(aes(DateTime, Temp/2), color="red")+
  scale_y_continuous(limits=c(5,11), sec.axis = sec_axis(~.*2, name=expression("Temp C")))



## Start building output spreadsheet
Output <- O2[,c("DateTime","Temp","DO")]
## Round DateTime to the nearest 5 minute interval to help with integration
Output$DateTime <- floor_date(Output$DateTime, "5 minutes")



#########
## CO2 ##
#########
CO2raw <- read.csv("2018_08_23_BlaineCreek_eosense.csv", header=T)
## Fix datetime
CO2raw$DateTime <- lubridate::mdy_hms(paste(CO2raw$Date, CO2raw$Time))
CO2raw$DateTime <- as.POSIXct(as.character(CO2raw$DateTime), format="%Y-%m-%d %H:%M:%S")
head(CO2raw)
## Fix names
CO2rawdata <- CO2raw[,c("Low.Range.CO2..ppm.","Temp..C.","Hi.Range.CO2..ppm.","DateTime")]
names(CO2rawdata) <- c("CO2ppm_LowR", "TempC","CO2ppm_HighR","DateTime")
CO2rawdata[,1:3] <- apply(CO2rawdata[,1:3],2,function(x) as.numeric(as.character(x)))
head(CO2rawdata)

sapply(CO2rawdata, class)

CO2 <- CO2rawdata

## Visualize
ggplot(CO2, aes(DateTime, CO2ppm_LowR)) +
  geom_point()+
  geom_point(aes(DateTime, TempC*50), color="red")+
  scale_y_continuous(limits=c(0,1500), sec.axis = sec_axis(~./50, name=expression("Temp C")))

## Calculate actual CO2 conc - Eosense
names(CO2)
## Import filled 5 minute Baro from FLBS North Shore, MT station
BP <- read.csv("Baro_1sec_BlaineCreek_2018_08_16.csv", header=T)
BP$DateTime <- as.POSIXct(as.character(BP$DateTime), format="%Y-%m-%d %H:%M:%S")

## Convert barometric pressure into atm from mmHg to atm
alt <- 901 #meter elevation at Blaine
BP$Baro_atm <- (BP$Baro_mmHg*exp( (-9.80665*0.0289644*alt)/(8.31447*(273.15+15)) ) )/760

CO2_BP <- merge(CO2, BP, by="DateTime")
names(CO2_BP)
head(CO2_BP)

## Rolling average of the previous 30 seconds
## Then select only time points at 5 minute intervals
CO2_BP <- CO2_BP[,c("DateTime","CO2ppm_LowR","TempC","Baro_atm")]

CO2_MA <- CO2_BP %>%
  tq_mutate(select = c(CO2ppm_LowR, TempC, Baro_atm), mutate_fun = rollapply,
            width = 8, align="right", FUN = mean, na.rm = TRUE,
            col_rename = c("CO2ppm_LowR_MA","TempC_MA","Baro_atm_MA"))

## Select only data from every 5 minutes (DO timeseries)
CO2_sub <- CO2_MA

CO2_sub$DateTime <- floor_date(CO2_sub$DateTime, unit = "5 seconds")
names(CO2_sub)
CO2_sub <- merge(Output[,c("DateTime","Temp")], CO2_sub[,c(1,5,6,7)], by="DateTime")
## If duplicates average across them
CO2_sub <- CO2_sub %>% group_by(as.numeric(CO2_sub$DateTime)) %>% summarise_all(funs(mean))
names(CO2_sub)
CO2_sub <- CO2_sub[,c(2,4:6)]
colnames(CO2_sub) <- c("DateTime","CO2ppm_LowR","TempC","Baro_atm")

CO2_sub$Caq_mmolm3 <- (CO2_sub$CO2ppm_LowR/38)*10^6






Converting2DissolvedConc_v2 <- function(x){
  #PV = nRT, gas conversion factor in L/mol
  #x$V <- (1*0.08206*(x$TempC+273.15))/x$Baro_atm
  
  #Hcc Henry's law dimensionless proportionality constant
  #x$kH_CO2 <- 29.41*exp(2400*((1/(x$TempC+273.15) - 1/(298.15)))) ## Van't Hoff equation
  #x$Hcc_CO2 <- 1/(x$kH_CO2)
  
  x$Caq_mmolm3 <- (x$CO2ppm_lowR/38)*10^6
  return(x)
}

## Now calculate CO2 concentrations
Converting2DissolvedConc <- function(x){
  #PV = nRT, gas conversion factor in L/mol
  x$V <- (1*0.08206*(x$TempC+273.15))/x$Baro_atm
  
  #Hcc Henry's law dimensionless proportionality constant
  x$kH_CO2 <- 29.41*exp(2400*((1/(x$TempC+273.15) - 1/(298.15)))) ## Van't Hoff equation
  x$Hcc_CO2 <- 1/(x$kH_CO2)
  
  #Low Range CO2
  #Converting Cg into umol/L (Cgas)
  x$Cgas_low <- x$CO2ppm_LowR/(x$V)
  #Finding concentration in aqueous soln (Caq) in mol/L
  x$Caq_molL_low <- x$Cgas_low/35
  x$Caq_mmolL_low <- x$Caq_molL_low*1000
  x$pCO2_uatm_Low <- x$Caq_molL_low*x$Baro_atm*1000
  #converting to mg/L
  #x$Caq_mgL_low <- x$Caq_molL_low*44.01*1000
    
  #high Range CO2
  #Converting Cg into umol/L (Cgas)
  #x$Cgas_high <- x$CO2ppm_HighR/(x$V)
  #Finding concentration in aqueous soln (Caq) in mol/L
  #x$Caq_molL_high <- x$Cgas_high*x$Hcc_CO2
  #x$pCO2_uatm_high <- x$Caq_molL_high*x$Baro_atm*1000
  #converting to mg/L
  #x$Caq_mgL_high <- x$Caq_molL_high*44.01*1000
  
  return(x)
}  

CO2_conc <- Converting2DissolvedConc_v2(CO2_sub)
names(CO2_conc)

## Select Date Range of interest
CO2_conc_diel <- subset(CO2_sub, DateTime >= "2018-08-22 05:30:00" & DateTime <= "2018-08-23 06:30:00")
names(CO2_conc_diel)

ggplot(CO2_conc_diel, aes(DateTime, (Caq_mmolm3/10^6)))+geom_point()

## Combine with previous Output
Output <- merge(Output, CO2_conc_diel, by="DateTime")


###############
## pH Sensor ##
###############
pHraw <- read.csv("pH_logfile_2018_08_24.csv", header=T)

pHraw$DateTime <- lubridate::mdy_hms(paste(as.character(pHraw$Date),as.character(pHraw$Time)))
pHraw$DateTime <- as.POSIXct(as.character(pHraw$DateTime), format="%Y-%m-%d %H:%M:%S",tz = "UTC")

## convert mV to pH using calibration values
buffer <- c(4.01, 7.00, 10.01)
mV <- c(1.17, 2.00, 2.83)
calib <- as.data.frame(t(rbind(buffer, mV)))
calib_sum <- summary(lm(calib$buffer ~ calib$mV))

pHraw$pH <- (calib_sum$coefficients[2,1]*pHraw$Voltage) + calib_sum$coefficients[1,1]

## Visualize
sapply(pHraw, class)

pH_5min <- pHraw %>%
  group_by(DateTime = cut(DateTime, breaks="5 min")) %>%
  dplyr::summarise(pH = mean(pH))
pH_5min$DateTime <- as.POSIXct(as.character(pH_5min$DateTime), format="%Y-%m-%d %H:%M:%S")
pH_5min$DateTime <- lubridate::round_date(pH_5min$DateTime, "5 minutes") 

ggplot(pH_5min, aes(DateTime, pH))+
  geom_point()


## Combine with previous Output
Output_ts5 <- merge(Output, pH_5min, by="DateTime", all.x = T)

###############
## Light
############
Light <- read.csv("NLDAS_5M_BlaineCreek.csv", header=T)
Light <- Light[,c("DateTime","PAR_filled")]
Light$DateTime <- as.POSIXct(as.character(Light$DateTime), format="%Y-%m-%d %H:%M", tz='UTC')
attributes(Light$DateTime)$tzone <- "America/Denver"


Output_ts5 <- merge(Output_ts5, Light, by="DateTime", all.x=T)

#write.csv(Output_ts5, "Compiled_5min_ts_Blaine.csv")

#######################################################
## Picarro Processed Isotope Data & Other Parameters ##
#######################################################
iso <- read.csv("2018_08_23_BlaineDiel_QAQC.csv", header=T)

iso <- iso[-which(iso$Syringe == "31"),]

## Separate SampleID
sapply(iso, class)
iso$SampleID <- as.character(iso$SampleID)

iso <- iso %>% 
  separate(SampleID, c("Site","Date","Time"))

iso$Date <- revalue(iso$Date, c("822" = "2018-08-22","823"="2018-08-23"))
iso$Time <- as.numeric(iso$Time)
iso$Time <- ifelse(iso$Time < 1000, yes= paste(substring(iso$Time, 1, 1), substring(iso$Time, 2, 3), sep = ":"),
                  no = ifelse(iso$Time >= 1000, yes=paste(substring(iso$Time, 1, 2), substring(iso$Time, 3, 4), sep = ":"), no=NA))
iso$Time <- revalue(iso$Time, c("0:"="0:00"))

iso$SampleDateTime <- lubridate::ymd_hm(paste(iso$Date,iso$Time))
iso$SampleDateTime <- as.POSIXct(as.character(iso$SampleDateTime), format="%Y-%m-%d %H:%M:%S",tz = "UTC")


## Calculate actual CO2 conc - Picarro
names(iso)

## For Picarro Data -- already has notes data attached
notes_iso <- iso
## Convert barometric pressure into atm from inHg to mmHg to atm
alt <- 901 #meter elevation at Blaine
notes_iso$Baro_atm <- ( notes_iso$Baro_inHg*25.4*exp( (-9.80665*0.0289644*alt)/(8.31447*(273.15+15)) ) )/760


CO2calc<-function(x) {
  x$kH_CO2<-29* exp( 2400 * ((1/(x$equiltemp+273.15)) - 1/(298.15)) ) #temp correction for Henry's law
  x$CO2_mol_water<- x$CO2_vol_air  / x$kH_CO2
  x$CO2_mol_air<- x$CO2_vol_air/ ((x$equiltemp+273.15)*0.08026/x$pressure)
  
  x$CO2_mmolm3<-(x$CO2_mol_water*x$volwater + x$CO2_mol_air*x$volair)/100 # µmol/L or mmol/m3
  x$CO2_uatm <- (x$CO2_mmolm3/1000)*x$pressure
  return(x)
  
}

## Need to set certain columns to desired values if all the same
## If each sample is different, simple set the column to a column specified in the metadata file
names(notes_iso)

notes_iso$CO2_vol_air <- notes_iso$CO2
notes_iso$equiltemp <- notes_iso$WaterTemp_C
notes_iso$volair <- notes_iso$air_mL
notes_iso$volwater <- notes_iso$H2O_mL
notes_iso$pressure <- notes_iso$Baro_atm

## Run function on datafile
notes_iso <- CO2calc(notes_iso)


## Summarize
notes_iso$DateTime <- notes_iso$SampleDateTime

iso_melt <- melt(notes_iso[,c("CO2","delCO2","CH4","delCH4","CO2_mmolm3","DateTime","Site")], id.vars=c("DateTime","Site"))
ggplot(iso_melt, aes(DateTime, value, color=Site))+
  geom_point()+
  facet_wrap(~variable, ncol=2, scales = "free")

iso_mean <- dcast(iso_melt, DateTime + Site ~ variable, fun.aggregate = mean)
iso_sd <- dcast(iso_melt, DateTime + Site ~ variable, fun.aggregate = sd)

iso_mean_melt <- melt(iso_mean[,c("CO2","delCO2","CH4","delCH4","CO2_mmolm3","DateTime","Site")], id.vars=c("DateTime","Site"))
iso_sd_melt <- melt(iso_sd[,c("CO2","delCO2","CH4","delCH4","CO2_mmolm3","DateTime","Site")], id.vars=c("DateTime","Site"))

iso_mean_sd <- merge(iso_mean_melt, iso_sd_melt, by=c("DateTime","Site","variable"))
colnames(iso_mean_sd)[4:5] <- c("Mean","SD")

sapply(iso_mean_sd, class)

## Visualize
ggplot(iso_mean_sd, aes(DateTime, Mean, color=Site))+
  geom_point(size=3)+
  geom_errorbar(aes(ymax=Mean-SD, ymin=Mean+SD), width=.2, position=position_dodge(.9))+
  facet_wrap(~variable, ncol=2, scales = "free")+
  theme(axis.title = element_blank(),
        strip.background = element_rect(fill = "white", color="black",linetype = "solid"),
        panel.border = element_rect(color = "black", linetype = "solid"),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))
  #scale_color_manual(values = c("B" = "red3"))
  #scale_x_datetime(breaks = iso_mean_sd$DateTime)

plot_grid(
  ggplot(iso_mean_sd[which(iso_mean_sd$variable == "CO2_mmolm3"),], aes(DateTime, Mean))+
    geom_point(size=4)+
    geom_errorbar(aes(ymax=Mean-SD, ymin=Mean+SD), width=0.2, position=position_dodge(.9))+
    xlab("")+
    ylab(expression('CO'[2]* ~ (mmol/m^{3})))+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_text(size=16), axis.text.y = element_text(size=14)),
  
  ggplot(iso_mean_sd[which(iso_mean_sd$variable == "delCO2"),], aes(DateTime, Mean))+
    geom_point(size=4)+
    geom_errorbar(aes(ymax=Mean-SD, ymin=Mean+SD), width=0.2, position=position_dodge(.9))+
    xlab("22-23 August 2018")+
    ylab(expression(delta^13*'C-CO'[2]))+
    theme(axis.title = element_text(size=16),
          axis.text = element_text(size=14))+
    scale_x_datetime(date_labels = "%H:%M"),
  ncol=1)






## Export
#write.csv(iso_mean_sd,"2018_08_23_Blaine_iso_CO2conc.csv")









#####################################################################################
## Visualize all together
######################################################################################

time_breaks <- c(as.POSIXct("2018-08-22 14:00:00"),
                 as.POSIXct("2018-08-22 18:00:00"),
                 as.POSIXct("2018-08-22 22:00:00"),
                 as.POSIXct("2018-08-23 02:00:00"),
                 as.POSIXct("2018-08-23 06:00:00"))

###############################
## pH sensor versus handheld ##
###############################

## merge to correlate
notes_minus5 <- notes
#notes_minus5$DateTime <- notes_minus5$DateTime - minutes(5)

pH_comparison <- merge(notes_minus5[,c("DateTime","pH")], pH_10min[,c("DateTime","pH")], by="DateTime")
names(pH_comparison) <- c("DateTime","Handheld_pH","Sensor_pH")

plot_grid(
  ggplot(pH_10min, aes(DateTime, pH))+
  geom_point(size=2,color="grey40")+
  geom_point(data=notes, aes(DateTime, pH), color="red2", size=3)+
  theme(axis.text.x = element_text(angle=45, hjust=1)),

ggplot(pH_comparison, aes(Handheld_pH, Sensor_pH))+
  geom_point(size=3)+
  scale_x_continuous(limits = c(7.8, 8.7))+
  scale_y_continuous(limits = c(7.8, 8.7))+
  geom_abline(slope=1, size=1, linetype="dashed")+
  xlab("Handheld pH")+
  ylab("Sensor pH"),

ncol=2)


#######################
## Gases and isotopes
######################



## CO2 and CO2 isotopic fractionation
ggplot(iso_mean_sd[which(iso_mean_sd$variable %in% c("CO2_conc","delCO2")),], aes(DateTime, Mean, color=Site))+
  geom_point(size=3, color="black")+
  geom_errorbar(aes(ymax=Mean-SD, ymin=Mean+SD), width=.2, position=position_dodge(.9), color="black")+
  facet_wrap(~variable, scales = "free", ncol=1)+
  theme(axis.title = element_blank(),
        strip.background = element_rect(fill = "white", color="black",linetype = "solid"),
        panel.border = element_rect(color = "black", linetype = "solid"),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))
  #ylab(expression(delta^13*'C-CO'[2]))

## CH4 and CH4 isotopic fractionation
ggplot(iso_mean_sd[which(iso_mean_sd$variable %in% c("CH4","delCH4")),], aes(DateTime, Mean, color=Site))+
  geom_point(size=3, color="black")+
  geom_errorbar(aes(ymax=Mean-SD, ymin=Mean+SD), width=.2, position=position_dodge(.9), color="black")+
  facet_wrap(~variable, scales = "free", ncol=1)+
  theme(axis.title = element_blank(),
        strip.background = element_rect(fill = "white", color="black",linetype = "solid"),
        panel.border = element_rect(color = "black", linetype = "solid"),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))
#ylab(expression(delta^13*'C-CO'[2]))


iso_plot <- merge(notes, iso_mean_sd[which(iso_mean_sd$variable %in% c("delCO2")),], by="DateTime", all.x = TRUE)
sapply(iso_plot,class)

plot_grid(
ggplot(iso_plot, aes(DateTime, Mean))+geom_point(size=3)+
    geom_errorbar(aes(ymax=Mean-SD, ymin=Mean+SD), width=.2, position=position_dodge(.9))+
    ylab("delCO2")+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
ggplot(iso_plot, aes(DateTime, pH))+geom_point(size=3)+
  geom_point(data=pH_10min, aes(DateTime,pH), color="grey20"),
  ncol=1)

## DO and pCO2 during diel
ggplot(O2, aes(DateTime, DO))+
  geom_point(color="blue")+
  ylab("Dissolved Oxygen (mg/L)")+
  geom_point(data=CO2_conc_diel, aes(DateTime,pCO2_uatm_Low/125))+
  scale_y_continuous(limits=c(3,11), sec.axis = sec_axis(~.*125, name=expression(paste("CO2 (",mu,"uatm)"))))


## DO and pCO2 full
ggplot(O2raw, aes(DateTime, DO))+
  geom_point(color="blue")+
  ylab("Dissolved Oxygen (mg/L)")+
  geom_point(data=CO2_conc, aes(DateTime,pCO2_uatm_Low/140))+
  scale_y_continuous(limits=c(3,11), sec.axis = sec_axis(~.*140, name=expression(paste("CO2 (",mu,"uatm)"))))

