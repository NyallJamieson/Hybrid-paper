# SDE-BP Coefficient of Variation

# Set initial parameter estimates
a <- 0.089
b <- 1.088
h <- 0.041
g <- a+b
G <- 62.792
r <- 32.072
theta <- (a^2+2*a*b+4*a*G*h-2*a*h+b^2-2*b*h+h^2)^0.5
A1 <- -0.5*(g+h+theta)
A2 <- -0.5*(g+h-theta)
B1 <- g-h+theta
B2 <- g-h-theta
B3 <- B1*B2

# Initialied vectors
times_vec_T1 <- c()
eps_vec_T1 <- c()

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

# Run experiment over ICS
for (n in 1:100){
  print(n)
  
  MEAN <- n*total_pop_mean
  VAR <- n*total_pop_var
  sd_over_mean_bp <- VAR^0.5/MEAN
  
  # Coefficient of Variation for BP
  
  dt <- 0.1
  t <- seq(0,150,by=dt)
  
  trajectory <- function(y){
    repeat {
      # Simulate feller simplified (allen only)
      # conditioned on not going extinct
      dt <- 0.1
      t <- seq(0,150,by=dt)
      l <- rep(0,length(t))
      m <- rep(0,length(t))
      l[1] <- 1
      noise1 <- rnorm(n=length(t),mean=0,sd=sqrt(dt))
      noise2 <- rnorm(n=length(t),mean=0,sd=sqrt(dt))
      
      for (i in 2:length(t)){
        # allen method
        A <- ((a+b)*l[i-1]+h*(G^2*(1+1/r)+G)*m[i-1])*(a*l[i-1]+h*m[i-1])-(a*l[i-1]+h*G*m[i-1])^2
        B <- sqrt(A)
        C <- (2*a+b)*l[i-1]+h*(G^2*(1+1/r)+G+1)*m[i-1]
        D <- sqrt(C+2*B)
        l[i] <- l[i-1] - dt*(a+b)*l[i-1] + dt*h*G*m[i-1] + ((a+b)*l[i-1]+h*(G^2*(1+1/r)+G)*m[i-1]+B)*noise1[i]/D - (a*l[i-1]+h*G*m[i-1])*noise2[i]/D
        m[i] <- m[i-1] + dt*a*l[i-1] - dt*h*m[i-1] - (a*l[i-1]+h*G*m[i-1])*noise1[i]/D + (a*l[i-1]+h*m[i-1]+B)*noise2[i]/D
        
        if (is.na(l[i])){
          l[i] <- 0}
        if (is.na(m[i])){
          m[i] <- 0}
        if (l[i]<0){
          l[i] <- 0}
        if (m[i]<0){
          m[i] <- 0}
        if (l[i]==0 && m[i]==0){
          break}
      }
      
      if (l[i]!=0 && m[i]!=0){
        return(data.frame("L"=l,"M"=m))
      }
    }
  }
  
  # Simulate 100000 runs using parallel::mclapply
  results <- parallel::mclapply(1:100000,trajectory, mc.cores = 5)
  
  # Separate L and M columns
  L_columns <- do.call(cbind, lapply(results, `[[`, "L"))
  M_columns <- do.call(cbind, lapply(results, `[[`, "M"))
  
  #
  means1 <- rowMeans(L_columns)
  means2 <- rowMeans(M_columns)
  var1 <- c()
  var2 <- c()
  covar <- c()
  
  for (i in 1:length(t)){
    vals1 <- as.numeric(L_columns[i,])
    vals2 <- as.numeric(M_columns[i,])
    var1 <- append(var1,var(vals1))
    var2 <- append(var2,var(vals2))
  }
  
  for (i in 1:length(times)){
    covar <- append(covar,cov(as.numeric(L_columns[i,]),as.numeric(M_columns[i,])))}
  
  means <- means1+means2
  sd <- sqrt(var1+var2+2*covar)
  frac <- sd/means
  
  F <- (sd_over_mean_bp/frac)[251:1501]
  m <- F[1251]
  eps <- seq(0.00001,0.1,by=0.00001)
  
  q <- F[1251]
  CTimes <- c()
  
  for (i in 1:length(eps)){
    CTimes <- append(CTimes,out1[1501-length(which(abs(F-q)<eps[i])),1])}
  grad <- (CTimes[length(CTimes)]-CTimes[1])/(eps[length(eps)]-eps[1])
  
  dis <- c()
  for (i in 1:length(eps)){
    dis <- append(dis,CTimes[1]+grad*eps[i]-CTimes[i])}
  t2 <- which(dis==max(dis))
  
  times_vec_T1 <- append(times_vec_T1,CTimes[t2])
  eps_vec_T1 <- append(eps_vec_T1,eps[t2])
  print(times_vec_T1)
}
  
  
# Save plot with 600 DPI resolution
tiff("VaryIC2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(seq(1,100,by=1),times_vec_T1,col="red",type="l",lwd=3,xlab="Deposited dose",ylab="Time in hours",main=expression("T"[2]~" estimate"))

# Close the device
dev.off()
