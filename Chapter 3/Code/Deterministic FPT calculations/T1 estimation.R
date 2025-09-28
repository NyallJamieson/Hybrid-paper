# sde simulation
a <- 0.089
b <- 1.088
h <- 0.041
m <- 62.792
G <- m
r <- 32.072
dt <- 0.1
t <- seq(0,150,by=dt)

simulate_extinction <- function(y){
  time <- 0
  extinct <- 0
  l <- rep(0, length(t))
  m <- rep(0, length(t))
  l[1] <- 1
  m[1] <- 0
  noise1 <- rnorm(n=length(t), mean=0, sd=sqrt(dt))
  noise2 <- rnorm(n=length(t), mean=0, sd=sqrt(dt))

  for (i in 2:length(t)){
    # Allen method
    A <- ((a+b)*l[i-1] + h*(G + (1/r+1)*G^2)*m[i-1]) * (a*l[i-1] + h*m[i-1]) - (a*l[i-1] + h*G*m[i-1])^2
    B <- sqrt(A)
    C <- (2*a+b)*l[i-1] + h*(G + (1/r+1)*G^2 + 1)*m[i-1]
    D <- sqrt(C + 2*B)

    l[i] <- l[i-1] - dt*(a+b)*l[i-1] + dt*h*G*m[i-1] + ((a+b)*l[i-1] + h*(G + (1/r+1)*G^2)*m[i-1] + B)*noise1[i]/D - (a*l[i-1] + h*G*m[i-1])*noise2[i]/D
    m[i] <- m[i-1] + dt*a*l[i-1] - dt*h*m[i-1] - (a*l[i-1] + h*G*m[i-1])*noise1[i]/D + (a*l[i-1] + h*m[i-1] + B)*noise2[i]/D

    if (is.na(l[i])) l[i] <- 0
    if (is.na(m[i])) m[i] <- 0
    if (l[i] < 0) l[i] <- 0
    if (m[i] < 0) m[i] <- 0

    if (l[i]==0 && m[i]==0){
      extinct <- 1
      time <- t[i]
      break
    }

  }

  return(time)
}

# Run simulations in parallel
n_simulations <- 100000  # Number of trajectories
n_cores <- 2  # Number of cores to use

extinction_times <- unlist(parallel::mclapply(
  X = 1:n_simulations,               # Simulations index (1 to 500)
  FUN = simulate_extinction,         # The function to run
  mc.cores = n_cores                 # Number of cores to use
))

extinction_times <- extinction_times[which(extinction_times>0)]

tally <- c()
for (i in t){
  tally <- append(tally,length(which(extinction_times<i)))
}

cdf_sde <- tally / n_simulations

# Compute PDF from CDF
pdf_sde <- diff(cdf_sde) / dt

# Adjust time for PDF (midpoints of t)
t_mid <- t[-1]  # Remove the first element to align with diff(prob)

# Now BP probability of extinction code
# Extinction G NB
model <- function(t,y,parms){
  dy1 <- -(a+b)*y[1]+a*y[2]+b
  dy2 <- -h*y[2]+h*(r/(m+r-y[1]*m))^(r)
  list(c(dy1,dy2))
}
yini <- c(y1=0,y2=0)
times <- seq(from=0,to=150,by=dt)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)

cdf_bp <- out1[,2]
###

# Compute PDF from CDF
pdf_bp <- diff(cdf_bp) / dt

diff_in_prob <- abs(pdf_sde - pdf_bp)

# Generate the sequence for eps
eps <- seq(0.00001, 0.1, by=0.00001)

# Function to process a single epsilon
compute_time <- function(epsilon) {
  # Find the first time point where diff_in_prob is less than the epsilon
  t_mid[which(diff_in_prob < epsilon)[1]] + dt
}

# Use mclapply for parallel processing
TIMES <- unlist(parallel::mclapply(
  X = eps,                # The list of epsilons
  FUN = compute_time,     # The function to apply
  mc.cores = 2  # Use all but one core
))

# Now produce plots
plot(eps,TIMES,col="red",lwd=2,type="l",xlab=expression(epsilon[1]),ylab=expression(T[1]),main="Threshold time for extinction condition")
grad <- (TIMES[length(TIMES)]-TIMES[1])/(eps[length(eps)]-eps[1])
dis <- c()
for (i in 1:length(eps)){
  dis <- append(dis,TIMES[1]+grad*eps[i]-TIMES[i])}
t <- which(dis==max(dis))
lines(c(eps[1],eps[length(eps)]),c(TIMES[1],TIMES[length(TIMES)]),type="l",lty=2)
lines(seq(0,2*eps[t],by=0.00000001),TIMES[t]-grad*(eps[t]-seq(0,2*eps[t],by=0.00000001)),lty=2)
symbols(eps[t],TIMES[t],circles=1,add=TRUE,inches=0.1)

eps[t]
TIMES[t]

# Save plot with 600 DPI resolution
tiff("Fig3.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(t_mid,pdf_bp,xlim=c(0,10),lwd=2,type="l",col="red",xlab="Time in hours",ylab="Probability of extinction",main="Probability of extinction over time")
lines(t_mid,pdf_sde,lwd=2,col="blue")

# Increase legend text size using cex
legend("topright", legend = c(expression(p[BP]), expression(p[SDE])), lwd = 2, col = c("red", "blue"), cex = 1.5)

# Close the device
dev.off()
