# Set parameter estimates
a <- 0.089
b <- 1.088
h <- 0.041
m <- 62.792
G <- m
r <- 32.072
dt <- 0.1
t <- seq(0,150,by=dt)

# Initialise vectors
times_vec_T1 <- c()
eps_vec_T1 <- c()

# Extinction G NB
model <- function(t,y,parms){
  dy1 <- -(a+b)*y[1]+a*y[2]+b
  dy2 <- -h*y[2]+h*(r/(m+r-y[1]*m))^(r)
  list(c(dy1,dy2))
}
yini <- c(y1=0,y2=0)
times <- seq(from=0,to=150,by=0.1)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)

# IC to vary over
initial_legionella <- seq(1,100,by=1)

# Now BP probability of extinction code

for (n in initial_legionella){
  print(n)
  cdf_bp <- out1[,2]^n
  # Compute PDF from CDF
  pdf_bp <- diff(cdf_bp) / dt
  
  simulate_extinction <- function(y){
    time <- 0
    extinct <- 0
    l <- rep(0, length(t))
    m <- rep(0, length(t))
    l[1] <- n
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
  n_cores <- 5  # Number of cores to use
  
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
  
  diff_in_prob <- abs(pdf_sde - pdf_bp)
  
  
  # Generate the sequence for eps
  eps <- seq(0.00001, 0.1, by=0.00001)
  
  compute_time <- function(epsilon){
    # Find the index of the peak
    peak_index <- which.max(diff_in_prob)
    # Subset the vector after the peak
    post_peak_values <- diff_in_prob[(peak_index + 1):length(diff_in_prob)]
    # Find the first index after the peak where the value falls below epsilon
    first_below_epsilon <- which(post_peak_values < epsilon)[1]
    # Calculate the corresponding time index in the original vector
    if (!is.na(first_below_epsilon)) {
      time_point <- peak_index + first_below_epsilon
      cat("The curve first drops below epsilon at time point:", time_point, "\n")
    } else {
      cat("The curve does not drop below epsilon after the peak.\n")
    }
    return(t_mid[time_point]+dt)
  }
  
  # Use mclapply for parallel processing
  TIMES <- unlist(parallel::mclapply(
    X = eps,                # The list of epsilons
    FUN = compute_time,     # The function to apply
    mc.cores = 5  # Use all but one core
  ))
  
  grad <- (TIMES[length(TIMES)]-TIMES[1])/(eps[length(eps)]-eps[1])
  dis <- c()
  for (i in 1:length(eps)){
    dis <- append(dis,TIMES[1]+grad*eps[i]-TIMES[i])}
  tt <- which(dis==max(dis))
  
  times_vec_T1 <- append(times_vec_T1, TIMES[tt])
  eps_vec_T1 <- append(eps_vec_T1, eps[tt])
}


# Save plot with 600 DPI resolution
tiff("VaryIC1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(initial_legionella,times_vec_T1,col="red",type="l",lwd=3,xlab="Deposited dose",ylab="Time in hours",main=expression("T"[1]~"estimate"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("VaryIC1eps.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 1)  # Adjust margins (optional for better space)

# Produce plot
plot(seq(1,100,by=1), eps_vec_T1, col="red", type="l", lwd=3,
     xlab="Deposited dose",
     ylab=expression(epsilon[3]^{symbol("*")}),
     main=expression(epsilon[3]^"*"~" estimate"))

# Close the device
dev.off()
