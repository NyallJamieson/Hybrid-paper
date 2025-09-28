load("ker_G (1).RData")

cdf_G <- cumsum(ker_G$y)/sum(ker_G$y)
sample_from_dens_G <- function(n){
  u <- runif(n)
  sampled_values <- approx(cdf_G,ker_G$x,xout=u)$y
  return(sampled_values)
}

results_deterministic <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(results_deterministic) <- c("G","r","eps","T4","L","M")

xstart <- c(time=0,L=1,M=0)

tally <- 0

while (tally < 1000){
  print(tally)
  
  sample_for_estimating <- na.omit(round(sample_from_dens_G(100000)))
  if (length(sample_for_estimating)<1){
    while (length(sample_for_estimating)<1){
      sample_for_estimating <- as.numeric(na.omit(round(sample_from_dens_G(100000))))
    }
  }
  
  lambda <- 0.04050593
  
  G_r <- fitdistrplus::fitdist(as.numeric(sample_for_estimating),"nbinom",method="mle")
  G_r_est <- MASS::mvrnorm(n=1,mu=c(coef(G_r)[1],coef(G_r)[2]),Sigma=matrix(c(vcov(G_r)[1],vcov(G_r)[2],vcov(G_r)[3],vcov(G_r)[4]),2,2))
  G <- G_r_est[2]
  r <- G_r_est[1]
  
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
  
  ####################################################################
  # Coefficient of Variation
  
  g <- a+b
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
  times2 <- seq(from=0,to=150,by=0.1)
  out1 <- deSolve::ode(times=times2,y=yini,func=model,parms=NULL)
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
  times2 <- seq(from=0,to=150,by=0.1)
  yini <- c(y1=0,y2=0,y3=0,y4=0,y5=0,y6=0,y7=0,y8=0)
  out2 <- deSolve::ode(times=times2,y=yini,func=model,parms=NULL)
  total_pop_var <- out2[,2]+out2[,3]+out2[,4]+out2[,5]
  TIMES <- c()
  
  n <- 1
  MEAN <- n*total_pop_mean
  VAR <- n*total_pop_var
  sd_over_mean <- VAR^0.5/MEAN
  
  
  frac_past_201 <- sd_over_mean[202:length(sd_over_mean)]
  #####
  times <- c()
  eps <- seq(0.00001,0.1,by=0.00001)
  for (i in 1:length(eps)){
    times <- append(times,out1[length(which(frac_past_201>sd_over_mean[1501]+eps[i])),1]+20)}
  
  grad <- (times[length(times)]-times[1])/(eps[length(eps)]-eps[1])
  dis <- c()
  for (i in 1:length(eps)){
    dis <- append(dis,times[1]+grad*eps[i]-times[i])}
  
  t <- which(dis==max(dis))
  
  E <- eps[t]
  T <- times[t]
  
  L2 <- out1[length(which(times2<=T)),2]
  M2 <- out1[length(which(times2<=T)),3]
  
  tryCatch({
    # Code that might throw an error
    results_deterministic <- dplyr::bind_rows(results_deterministic,data.frame(G=G,r=r,eps=E,T4=T,L=L2,M=M2))
    
    # If no error occurs, increment the success counter
    tally <- tally  + 1
  }, error = function(e) {
    # Skip this iteration on error
    message("Skipping due to error: ", e$message)
  })
  
}

write.csv(results_deterministic,"sens_T4_deterministic_G_r.csv")

results_stochastic <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(results_stochastic) <- c("G","r","time")

stoch_sample_size <- 1000

stoch <- function(y){
  
  G <- results_deterministic[y,1]
  r <- results_deterministic[y,2]
  L2 <- results_deterministic[y,5]
  M2 <- results_deterministic[y,6]
  
  frame <- data.frame("G" = rep(G, stoch_sample_size), "r" = rep(r, stoch_sample_size), "time" = NA)
  
  
  #################################################
  # NB Jump SDE
  
  # tau leap simulation for Jump process
  times3 <- c()
  
  while (length(times3)<stoch_sample_size){
    
    dt <- 0.1
    t <- seq(0,300,by=dt)
    l <- rep(0,length(t))
    m <- rep(0,length(t))
    l[1] <- 1
    
    for (i in 2:3001){
      
      e1 <- rpois(n=1,lambda=dt*b*l[i-1])
      e2 <- rpois(n=1,lambda=dt*a*l[i-1])
      e3 <- rpois(n=1,lambda=dt*h*m[i-1])
      if (e3==0){
        e4 <- 0}
      else {
        e4 <- sum(rnbinom(n=e3,mu=G,size=r))}
      
      l[i] <- l[i-1] - (e1+e2) + e4
      m[i] <- m[i-1] + e2 - e3
      
      if (l[i]<0){
        l[i] <- 0}
      if (m[i]<0){
        m[i] <- 0}
      
      if (l[i]>L2 && m[i]>M2){
        times3 <- append(times3,t[i])
        print(length(times3))
        break}
      
      
    }
  }
  frame[,3] <- times3
  return(frame)
}

results_stochastic <- dplyr::bind_rows(results_stochastic,as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:1000,stoch,mc.cores=45)),na.action="omit")))

saveRDS(results_stochastic,"sens_T4_stochastic_G_r.rds")
