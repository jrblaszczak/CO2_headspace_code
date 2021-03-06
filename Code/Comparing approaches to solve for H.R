library(seacarb)


####Fake data
AT=2000
co2=400
K1=0.000045
K2=0.000038
Kw=0.000004

####These equations are derived from Alk= Bicarbonate+2*Carbonate. They assume that H and OH are zero. 
####Only the polyroot equation can solve for 3rd order (or greater) polynomials. 

#Polyroot solution for solving H+ (signs as Kos et al.)
h_all <- polyroot(c((-2*K1*K2*co2),-(co2*K1),AT))####the final factor should be -1 by MTT's math
real<-Re(h_all)
h <-real[which(real>0)]
#h=3.103771e-05

#Polyroot solution for solving H+ (signs opposite Kos et al.)
h_all <- polyroot(c((2*K1*K2*co2),(co2*K1),-AT))
real<-Re(h_all)
h <-real[which(real>0)]
#h=3.103771e-05

#Matlab solution for solving H+
h_matlab<-(4*co2*K1*K2)/((co2*K1*(8*AT*K2+co2*K1))^0.5-co2*K1)
#h_matlab=3.103771e-05

#Solving by hand using Quadratic equation (Joanna)
h_quad=((-K1*co2)-((K1*co2)^2-4*(-AT*2*K1*K2*co2))^0.5)/(2*-AT)
#h_jo=3.103771e-05

#####The two options below are not working yet.

#Minimization function
f<-function(x, K1, K2, co2, AT) {
  -(2*K1*K2*co2)/x^2+(co2*K1)/x-AT
}
x<-0.00000012
xmin<-nlm(f, x,AT=2000,co2=400,K1=0.000045,K2=0.000038, print.level=2, steptol=0.0000001
)
xmin


#SeaCarb
sc<-carb(flag=4, co2, AT, S=20, kf="dg",k1k2="m06")
10^-(sc$pH)


