load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_lambda.RData")

cdf_lambda <- cumsum(ker_lambda$y)/sum(ker_lambda$y)
sample_from_dens_lambda <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_lambda,ker_lambda$x,xout=u)$y
  return(sampled_values)
}

results_deterministic <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(results_deterministic) <- c("lambda","eps","T1","L","M")

xstart <- c(time=0,L=1,M=0)

times_vec_T1 <- c()

tally <- 0

while (tally < 1000){
  
  error_check <- 0
  
  print(tally)
  
  phi <- 0.249
  lambda=sample_from_dens_lambda(1)
  A <- 177.95295
  w <- 0.192
  G <- round(A/(1+(A-1)*exp(-w*(1/lambda-1))))
  time <- c(2,24,48,72)
  logdose <- log(c(13584.15,588828.58,11207102.05,82907236.06))
  L0 <- 10^5
  yourdata2 <- data.frame(time,logdose)
  ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
  if (min(coef(ab_model))<0){
    while (min(coef(ab_model))<0){
      ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
    }
  }
  alpha <- unname(coef(ab_model)[1]*coef(ab_model)[2])
  beta <- unname(coef(ab_model)[1]*(1-coef(ab_model)[2]))
  
  a <- alpha
  b <- beta
  h <- lambda
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
  n_simulations <- 50000  # Number of trajectories
  n_cores <- 20  # Number of cores to use
  
  extinction_times <- unlist(parallel::mclapply(
    X = 1:n_simulations,               # Simulations index (1 to 500)
    FUN = simulate_extinction,         # The function to run
    mc.cores = n_cores                 # Number of cores to use
  ))
  
  extinction_times <- extinction_times[which(extinction_times>0)]
  
  total <- c()
  for (i in t){
    total <- append(total,length(which(extinction_times<i)))
  }
  
  cdf_sde <- total / n_simulations
  
  # Compute PDF from CDF
  pdf_sde <- diff(cdf_sde) / dt
  
  # Adjust time for PDF (midpoints of t)
  t_mid <- t[-1]  # Remove the first element to align with diff(prob)
  
  
  
  
  
  # now BP probability of extinction code
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
    mc.cores = 20  # Use all but one core
  ))
  
  # Remove NA values
  TIMES <- TIMES[!is.na(TIMES)]
  
  # Check if TIMES is empty or contains invalid values
  if (!is.numeric(TIMES) || length(TIMES[!is.na(TIMES)]) < 2) {
    message("Skipping iteration due to invalid TIMES.")
    next  # Skip this iteration and continue the while loop
  }
  
  grad <- (TIMES[length(TIMES)]-TIMES[1])/(eps[length(eps)]-eps[1])
  dis <- c()
  for (i in 1:length(eps)){
    dis <- append(dis,TIMES[1]+grad*eps[i]-TIMES[i])}
  tt <- which(dis==max(dis))
  
  E <- eps[tt]
  T <- TIMES[tt]
  
  # Mean
  model <- function(t,y,parms){
    dy1 <- -(a+b)*y[1]+h*G*y[2]
    dy2 <- a*y[1]-h*y[2]
    list(c(dy1,dy2))
  }
  yini <- c(y1=1,y2=0)
  times0 <- seq(from=0,to=150,by=0.1)
  out0 <- deSolve::ode(times=times0,y=yini,func=model,parms=NULL)
  
  tryCatch({
    # Code that might throw an error
    L2 <- out0[length(which(t<=T)),2]
    M2 <- out0[length(which(t<=T)),3]
    
  }, error = function(e) {
    # Skip this iteration on error
    error_check <<- 1
    message("dis error: ", e$message)
  })
  
  results_deterministic <- dplyr::bind_rows(results_deterministic,data.frame(lambda=lambda,eps=E,T1=T,L=L2,M=M2))
  tally <- tally + 1
}

write.csv(results_deterministic,"~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T1_lambda_deterministic.csv")
