library(seacarb)

AT=2000
co2=400
K1=0.000045
K2=0.000038
Kw=0.000004

#Polyroot solution for solving H+ (signs as Kos et al.)
h_all <- polyroot(c(-1*(2*K1*K2*co2),-1*(co2*K1),AT))
real<-Re(h_all)
h <-real[which(real>0)]
#h=3.103771e-05

#Polyroot solution for solving H+ (signs opposite Kos et al.)
h_all <- polyroot(c((2*K1*K2*co2),(co2*K1),-1*AT))
real<-Re(h_all)
h <-real[which(real>0)]
#h=3.103771e-05

#Matlab
h_matlab<-(4*co2*K1*K2)/((co2*K1*(8*AT*K2+co2*K1))^0.5-co2*K1)
#h_matlab=3.103771e-05

#Quadratic equation (Joanna)
h_jo=((-K1*co2)-((K1*co2)^2-4*(-AT*2*K1*K2*co2))^0.5)/(2*-AT)
#h_jo=3.103771e-05