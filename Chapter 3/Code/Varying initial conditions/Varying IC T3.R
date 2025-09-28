# Extinction G NB
model <- function(t,y,parms){
  a <- 0.089
  b <- 1.088
  h <- 0.041
  m <- 62.792
  G <- m
  r <- 32.072
  dy1 <- -(a+b)*y[1]+a*y[2]+b
  dy2 <- -h*y[2]+h*(r/(m+r-y[1]*m))^(r)
  list(c(dy1,dy2))
}
yini <- c(y1=0,y2=0)
times <- seq(from=0,to=150,by=0.1)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)

# Save plot with 600 DPI resolution
tiff("VaryIC3Example.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(out1[,1],out1[,2],xlim=c(0,20),ylim=c(0,1),xlab="Time in hours",ylab="Probability",main="Extinction probability",type="l",lwd=2)
ics <- c(1,2,3,4,5)
for (i in 1:length(ics)){
  lines(out1[,1],out1[,2]^ics[i],col=i+1,lwd=2)}

# Add legend
legend("bottomright",col=c(2,3,4,5,6),legend=c("L(0)=1","L(0)=2","L(0)=3","L(0)=4","L(0)=5"),lwd=2,cex=1.5)

# Close the device
dev.off()

# Start varying IC calculations now
TIMES <- c()
eps_values <- c()

for (n in 1:100){
  print(n)
  eps <- seq(0.00001,0.1,by=0.00001)
  OUT1 <- out1[,2]^n
  times <- c()
  
  for (i in 1:length(eps)){
    times <- append(times,out1[length(which(OUT1<OUT1[1501]-eps[i])),1])}
  
  eps <- eps[1:length(times)]
  grad <- (times[length(times)]-times[1])/(eps[length(eps)]-eps[1])
  dis <- c()
  
  for (i in 1:length(eps)){
    dis <- append(dis,times[1]+grad*eps[i]-times[i])}
  
  TIMES <- append(TIMES,times[which(dis==max(dis))])
  eps_values <- append(eps_values, eps[which(dis==max(dis))])
  }

# Save plot with 600 DPI resolution
tiff("VaryIC3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(seq(1,100,by=1),TIMES,col="red",type="l",lwd=3,xlab="Deposited dose",ylab="Time in hours",main=expression("T"[3]~" estimate"))

# Close the device
dev.off()

# Save plot with 600 DPI resolution
tiff("VaryIC3eps.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 1)  # Adjust margins (optional for better space)

# Produce plot
plot(seq(1,100,by=1), eps_values, col="red", type="l", lwd=3,
     xlab="Deposited dose",
     ylab=expression(epsilon[3]^{symbol("*")}),
     main=expression(epsilon[3]^"*"~" estimate"))

# Close the device
dev.off()