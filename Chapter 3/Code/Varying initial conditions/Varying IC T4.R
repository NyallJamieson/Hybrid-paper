# Set parameter estimates
a <- 0.089
b <- 1.088
h <- 0.041
m <- 62.792
G <- m
r <- 32.072
g <- a + b
theta <- (a^2+2*a*b+4*a*G*h-2*a*h+b^2-2*b*h+h^2)^0.5
A1 <- -0.5*(g+h+theta)
A2 <- -0.5*(g+h-theta)
B1 <- g-h+theta
B2 <- g-h-theta
B3 <- B1*B2
 
# Mean
model <- function(t,y,parms){
dy1 <- -(a+b)*y[1]+h*G*y[2]
dy2 <- a*y[1]-h*y[2]
list(c(dy1,dy2))
}
yini <- c(y1=1,y2=0)
times <- seq(from=0,to=150,by=0.1)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)
total_pop_mean <- out1[,2]+out1[,3]
 
# Variance
model <- function(t,y,parms){
dy1 <- -2*g*y[1]+h*G*(y[2]+y[3]) + 0.5*g*(B1*exp(t*A1)+B2*exp(t*A2))/theta + h*(G+(1+1/r)*G^2)*a*(exp(t*A2)-exp(t*A1))/theta
dy2 <- a*y[1]-(g+h)*y[2]+h*G*y[4] - 0.5*a*(B1*exp(t*A1)+B2*exp(t*A2))/theta - h*G*a*(exp(t*A2)-exp(t*A1))/theta
dy3 <- a*y[1]-(g+h)*y[3]+h*G*y[4] - 0.5*a*(B1*exp(t*A1)+B2*exp(t*A2))/theta - h*G*a*(exp(t*A2)-exp(t*A1))/theta
dy4 <- a*(y[2]+y[3])-2*h*y[4] + 0.5*a*(B1*exp(t*A1)+B2*exp(t*A2))/theta + h*a*(exp(t*A2)-exp(t*A1))/theta
dy5 <- -2*g*y[5]+h*G*(y[6]+y[7]) + g*B3*(exp(t*A1)+exp(t*A2))/(4*a*theta) + h*(G+(1+1/r)*G^2)*(B1*exp(t*A2)-B2*exp(t*A1))/(2*theta)
dy6 <- a*y[5]-(g+h)*y[6]+h*G*y[8] - B3*(exp(t*A1)+exp(t*A2))/(4*theta) - h*G*(B1*exp(t*A2)-B2*exp(t*A1))/(2*theta)
dy7 <- a*y[5]-(g+h)*y[7]+h*G*y[8] - B3*(exp(t*A1)+exp(t*A2))/(4*theta) - h*G*(B1*exp(t*A2)-B2*exp(t*A1))/(2*theta)
dy8 <- a*(y[6]+y[7])-2*h*y[8] + B3*(exp(t*A1)+exp(t*A2))/(4*theta) + h*(B1*exp(t*A2)-B2*exp(t*A1))/(2*theta)
list(c(dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8))
}
times <- seq(from=0,to=150,by=0.1)
yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0,y6=0,y7=0,y8=0)
out2 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)
total_pop_var <- out2[,2]+out2[,3]+out2[,4]+out2[,5]

# Initialise vectors
TIMES <- c()
eps_values <- c()

for (n in 1:100){
  print(n)
 
  MEAN <- n*total_pop_mean
  VAR <- n*total_pop_var
  sd_over_mean <- VAR^0.5/MEAN
 
  frac_past_201 <- sd_over_mean[202:length(sd_over_mean)]

  times <- c()
  eps <- seq(0.00001,0.1,by=0.00001)
  for (i in 1:length(eps)){
    times <- append(times,out1[length(which(frac_past_201>sd_over_mean[1501]+eps[i])),1]+20)
    }
 
  grad <- (times[length(times)]-times[1])/(eps[length(eps)]-eps[1])
  dis <- c()
  for (i in 1:length(eps)){
    dis <- append(dis,times[1]+grad*eps[i]-times[i])
    }
 
  TIMES <- append(TIMES,times[which(dis==max(na.omit(dis)))])
  eps_values <- append(eps_values,eps[which(dis==max(na.omit(dis)))])
  }
 
# Save plot with 600 DPI resolution
tiff("VaryIC4.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(seq(1,100,by=1),TIMES,type="l",lwd=2,col="red",xlab="Deposited dose",ylab="Time in hours",main=expression("T"[4]~" estimate"))

# Close the device
dev.off()

# Save plot with 600 DPI resolution
tiff("VaryIC4eps.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 1)  # Adjust margins (optional for better space)

# Produce plot
plot(seq(1,100,by=1), eps_values, col="red", type="l", lwd=3,
     xlab="Deposited dose",
     ylab=expression(epsilon[4]^{symbol("*")}),
     main=expression(epsilon[4]^"*"~" estimate"))

# Close the device
dev.off()