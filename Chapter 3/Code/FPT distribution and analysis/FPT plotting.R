# Read in the simulated FPT data
T1_gill <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T1_gill.csv")[,2]
T1_tau <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T1_tau.csv")[,2]
T2_gill <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T2_gill.csv")[,2]
T2_tau <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T2_tau.csv")[,2]
T2_sde <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T2_sde.csv")[,2]

# Obtain kernel densities for these datasets
dens_T1g <- density(T1_gill)
dens_T1t <- density(T1_tau)
dens_T2g <- density(T2_gill)
dens_T2t <- density(T2_tau)
dens_T2s <- density(T2_sde)

# Save plot with 600 DPI resolution
tiff("Fig4.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce FPT1 plot
hist(T1_gill,freq=FALSE,breaks=100,xlim=c(0,200),ylim=c(0,0.025),xlab="Time in hours",main=expression(T[SDE]^"*" * " distribution"))
lines(dens_T1t,col="blue",lwd=2)
lines(dens_T1g,col="red",lwd=2)

# Increase legend text size using cex
legend("topright",col=c("red","blue"),lwd=2,legend=c("Gillespie","Jump SDE"), cex = 1.5)

# Close the device
dev.off()

# Save plot with 600 DPI resolution
tiff("Fig4.2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce FPT2 plot
hist(T2_gill,freq=FALSE,breaks=150,xlim=c(0,200),ylim=c(0,0.025),xlab="Time in hours",main=expression(T[ODE]^"*" * " distribution"))
lines(dens_T2t,col="blue",lwd=2)
lines(dens_T2s,col="green",lwd=2)
lines(dens_T2g,col="red",lwd=2)

# Increase legend text size using cex
legend("topright",col=c("red","blue","green"),lwd=2,legend=c("Gillespie","Jump SDE","Gaussian SDE"), cex = 1.5)

# Close the device
dev.off()

# Low-dose incubation-period plot
a <- 0.089
b <- 1.088
h <- 0.041
m <- 62.792
G <- m
r <- 32.072

# Mean
model <- function(t,y,parms){
  dy1 <- -(a+b)*y[1]+h*G*y[2]
  dy2 <- a*y[1]-h*y[2]
  list(c(dy1,dy2))
}
yini <- c(y1=1,y2=0)
times <- seq(from=0,to=150,by=0.1)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)
T2 <- 49.9
T3 <- out1[length(which(out1[,2]<50661)),1]
T2to3 <- T3-T2

# BP from t=0 to T_{2}^{*}, then ODE to T_{3}^{*}
T3_BPtoT2_ODEtoT3 <- T2_gill + T2to3
dens_T3_BPtoT2_ODEtoT3 <- density(T3_BPtoT2_ODEtoT3)

# BP from t=0 to T_{1}^{*}, then SDE to T_{2}^{*}, then ODE to T_{3}^{*}
T3_BPtoT1_SDEtoT2_ODEtoT3 <- T2_sde + T2to3
dens_T3_BPtoT1_SDEtoT2_ODEtoT3 <- density(T3_BPtoT1_SDEtoT2_ODEtoT3)

# BP from t=0 to T_{1}^{*}, then SDE to T_{3}^{*}
T3_BPtoT1_SDEtoT3 <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T3_sde.csv")[,2]
dens_T3_BPtoT1_SDEtoT3 <- density(T3_BPtoT1_SDEtoT3)

# BP from t=0 to T_{3}^{*}
T3_gill <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T3_gill.csv")[,2]
dens_T3_gill <- density(T3_gill)

# Save plot with 600 DPI resolution
tiff("Fig5.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
hist(T3_gill,xlab="Time in hours",main="Low-dose incubation-period distribution",freq=FALSE,breaks=100,xlim=c(50,200))
lines(dens_T3_BPtoT2_ODEtoT3,col="green",lwd=2)
lines(dens_T3_BPtoT1_SDEtoT2_ODEtoT3,col="orange",lwd=2)
lines(dens_T3_BPtoT1_SDEtoT3,col="blue",lwd=2)
lines(dens_T3_gill,col="red",lwd=2)

# Increase legend text size using cex
legend("topright",col=c("red","blue","green","orange"), lwd=2,legend=c("Model 1","Model 2","Model 3","Model 4"), cex = 1.5)

# Close the device
dev.off()
