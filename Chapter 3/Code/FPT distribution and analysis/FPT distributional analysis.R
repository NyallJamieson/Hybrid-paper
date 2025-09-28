# models that I want to edit later and to source from a Chapter 1 script

# cumulative density function (c.d.f) of type III Burr distribution
pburr3 <- function(x,c,k,s){(1+(x/s)^(-c))^(-k)}
# probability density function (p.d.f) of type III Burr distribution
dburr3 <- function(x,c,k,s){(c*k/s)*(x/s)^(-c-1)*(1+(x/s)^(-c))^(-k-1)}
# c.d.f of type X Burr distribution
pburr10 <- function(x,c,k){(1-exp(-(x/c)^(2)))^(k)}
# p.d.f of type X Burr distribution
dburr10 <- function(x,c,k){2*k*x*c^(-2)*exp(-(x/c)^(2))*(1-exp(-(x/c)^(2)))^(k-1)}
# c.d.f of type XII Burr distribution
pburr12 <- function(x,c,k,s){1-(1+(x/s)^(c))^(1-k)}
# p.d.f of type XII Burr distribution
dburr12 <- function(x,c,k,s){c*(k-1)*s^(-1)*(1+(x/s)^(c))^(-k)*(x/s)^(c-1)}
# c.d.f of our derived Burr distribution
dburrNEW <- function(x,a,b,T){(x/b+a)*(T/x)^(a)*exp((T-x)/b)/(x*(1+(T/x)^(a)*exp((T-x)/b))^2)}
# p.d.f of our derived Burr distribution
pburrNEW <- function(x,a,b,T){(1+(T/x)^(a)*exp((T-x)/b))^(-1)}

FPT1 <- read.csv("T1_gill.csv")[,2]
FPT2 <- read.csv("T2_gill.csv")[,2]

# FPT 1 analysis

m1_1 <- mean(FPT1)
m2_1 <- var(FPT1)
m3_1 <- PerformanceAnalytics::skewness(FPT1)
m4_1 <- PerformanceAnalytics::kurtosis(FPT1,method="moment")

FPT1_model1 <- fitdistrplus::fitdist(FPT1,"gamma",method="mle")
FPT1_model2 <- fitdistrplus::fitdist(FPT1,"weibull",method="mle")
FPT1_model3 <- fitdistrplus::fitdist(FPT1,"lnorm",method="mle")
FPT1_model4 <- fitdistrplus::fitdist(FPT1,"burr3",method="mle",start=list(c=1,k=1,s=1))
FPT1_model5 <- fitdistrplus::fitdist(FPT1,"burr10",method="mle",start=list(c=50,k=0.5))
FPT1_model6 <- fitdistrplus::fitdist(FPT1,"burr12",method="mle",start=list(c=1,k=2,s=100))
FPT1_model7 <- fitdistrplus::fitdist(FPT1,"burrNEW",method="mle",start=list(a=1,b=1,T=1))

AIC(FPT1_model1)
AIC(FPT1_model2)
AIC(FPT1_model3)
AIC(FPT1_model4)
AIC(FPT1_model5)
AIC(FPT1_model6)
AIC(FPT1_model7)

which(c(AIC(FPT1_model1),AIC(FPT1_model2),AIC(FPT1_model3),AIC(FPT1_model4),AIC(FPT1_model5),AIC(FPT1_model6),AIC(FPT1_model7))==min(c(AIC(FPT1_model1),AIC(FPT1_model2),AIC(FPT1_model3),AIC(FPT1_model4),AIC(FPT1_model5),AIC(FPT1_model6),AIC(FPT1_model7))))
FPT1_model3

# FPT 1 analysis

m1_2 <- mean(FPT2)
m2_2 <- var(FPT2)
m3_2 <- PerformanceAnalytics::skewness(FPT2)
m4_2 <- PerformanceAnalytics::kurtosis(FPT2,method="moment")

FPT2_model1 <- fitdistrplus::fitdist(FPT2,"gamma",method="mle")
FPT2_model2 <- fitdistrplus::fitdist(FPT2,"weibull",method="mle")
FPT2_model3 <- fitdistrplus::fitdist(FPT2,"lnorm",method="mle")
FPT2_model4 <- fitdistrplus::fitdist(FPT2,"burr3",method="mle",start=list(c=10,k=1,s=1))
FPT2_model5 <- fitdistrplus::fitdist(FPT2,"burr10",method="mle",start=list(c=50,k=0.5))
FPT2_model6 <- fitdistrplus::fitdist(FPT2,"burr12",method="mle",start=list(c=1,k=2,s=100))
FPT2_model7 <- fitdistrplus::fitdist(FPT2,"burrNEW",method="mle",start=list(a=10,b=1,T=50))

AIC(FPT2_model1)
AIC(FPT2_model2)
AIC(FPT2_model3)
AIC(FPT2_model4)
AIC(FPT2_model5)
AIC(FPT2_model6)
AIC(FPT2_model7)

which(c(AIC(FPT2_model1),AIC(FPT2_model2),AIC(FPT2_model3),AIC(FPT2_model4),AIC(FPT2_model5),AIC(FPT2_model6),AIC(FPT2_model7))==min(c(AIC(FPT2_model1),AIC(FPT2_model2),AIC(FPT2_model3),AIC(FPT2_model4),AIC(FPT2_model5),AIC(FPT2_model6),AIC(FPT2_model7))))
FPT2_model3
