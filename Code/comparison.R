




Hall.CO2.simple<-(StmpCO2)
Kos.CO2.simple<-pCO2$`pCO2 simple headspace (micro-atm)`
plot(Hall.CO2.simple~Kos.CO2.complete)
model<-lm(Hall.CO2.simple~Kos.CO2.simple)
summary(model)




Hall.CO2.complete<-Carb2$CO2
Kos.CO2.complete<-pCO2$`pCO2 complete headspace (ppmv)`
plot(Hall.CO2.complete~Kos.CO2.complete)
model<-lm(Hall.CO2.complete~Kos.CO2.complete)
summary(model)


#Kos K1 and K2

Kos.K1=10^-(-126.34048+6320.813/(temp_equil+273.15)+19.568224*log(temp_equil+273.15))
Kos.K2=10^-(-90.18333+5143.692/(temp_equil+273.15)+14.613358*log(temp_equil+273.15))

K1calc<- function(temp_equil) { 10^( (-3404.71/(273.15+temp_equil)) + 14.844 -0.033*(temp_equil+273.15) )}
K2calc<- function(temp_equil) { 10^( (-2902.39/(273.15+temp_equil)) + 6.498 -0.0238*(temp_equil+273.15) )}

Hall.K1<-K1calc(temp_equil)
Hall.K2<-K2calc(temp_equil)

plot(Hall.K1~Kos.K2)
model<-lm(Hall.K1~Kos.K2)
summary(model)


