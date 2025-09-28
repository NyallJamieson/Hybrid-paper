load("~/PhD-work/Chapter 2/Data/Parameter distributions/ker_lambda.RData")

cdf_lambda <- cumsum(ker_lambda$y)/sum(ker_lambda$y)
sample_from_dens_lambda <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_lambda,ker_lambda$x,xout=u)$y
  return(sampled_values)
}

results_deterministic <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(results_deterministic) <- c("lambda","eps","T3","L","M")

xstart <- c(time=0,L=1,M=0)

tally <- 0

while (tally < 1000){
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
  
  # Extinction G NB
  model <- function(t,y,parms){
    dy1 <- -(a+b)*y[1]+a*y[2]+b
    dy2 <- -h*y[2]+h*(r/(m+r-y[1]*m))^(r)
    list(c(dy1,dy2))
  }
  yini <- c(y1=0,y2=0)
  times1 <- seq(from=0,to=150,by=0.1)
  out1 <- deSolve::ode(times=times1,y=yini,func=model,parms=NULL)
  
  ###
  
  eps <- seq(0.00001,0.1,by=0.00001)
  TIMES <- c()
  
  for (i in 1:length(eps)){
    TIMES <- append(TIMES,out1[length(which(out1[,2]<out1[1501,2]-eps[i])),1])}
  
  grad <- (TIMES[length(TIMES)]-TIMES[1])/(eps[length(eps)]-eps[1])
  dis <- c()
  for (i in 1:length(eps)){
    dis <- append(dis,TIMES[1]+grad*eps[i]-TIMES[i])}
  t <- which(dis==max(dis))
  
  E <- eps[t]
  T <- TIMES[t]
  
  # Mean
  model <- function(t,y,parms){
    dy1 <- -(a+b)*y[1]+h*G*y[2]
    dy2 <- a*y[1]-h*y[2]
    list(c(dy1,dy2))
  }
  yini <- c(y1=1,y2=0)
  times2 <- seq(from=0,to=150,by=0.1)
  out1 <- deSolve::ode(times=times2,y=yini,func=model,parms=NULL)
  total_pop_mean <- out1[,2]+out1[,3]
  
  L2 <- out1[length(which(times2<=T)),2]
  M2 <- out1[length(which(times2<=T)),3]
  
  tryCatch({
    # Code that might throw an error
    results_deterministic <- dplyr::bind_rows(results_deterministic,data.frame(lambda=lambda,eps=E,T3=T,L=L2,M=M2))
    
    # If no error occurs, increment the success counter
    tally <- tally  + 1
  }, error = function(e) {
    # Skip this iteration on error
    message("Skipping due to error: ", e$message)
  })
}

write.csv(results_deterministic, "~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T3_lambda_deterministic.csv")