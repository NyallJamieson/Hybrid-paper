results_deterministic <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(results_deterministic) <- c("alpha","beta","eps","T3","L","M")

xstart <- c(time=0,L=1,M=0)

tally <- 0

while (tally < 1000){
  print(tally)
  
  lambda <- 0.04050593
  A <- 177.95295
  w <- 0.192
  G <- 62
  time <- c(2,24,48,72)
  logdose <- log(c(13584.15,588828.58,11207102.05,82907236.06))
  L0 <- 10^5
  yourdata2 <- data.frame(time,logdose)
  ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
  gamma_pi <- MASS::mvrnorm(n=1,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2))
  if (min(gamma_pi)<0){
    while (min(gamma_pi)<0){
      ab_model <- minpack.lm::nlsLM(logdose~log(L0)-0.5*log((lambda-gamma)^2+4*lambda*Pi*gamma*G)+log(0.5*(lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma+((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))-0.5*(lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5)*exp(0.5*time*(-lambda-gamma-((lambda-gamma)^2+4*lambda*gamma*Pi*G)^0.5))),start=list(gamma=5,Pi=0.5),lower=c(0,0))
      gamma_pi <- MASS::mvrnorm(n=1,mu=c(coef(ab_model)[1],coef(ab_model)[2]),Sigma=matrix(c(vcov(ab_model)[1],vcov(ab_model)[2],vcov(ab_model)[3],vcov(ab_model)[4]),2,2))
    }
  }
  alpha <- unname(gamma_pi[1]*gamma_pi[2])
  beta <- unname(gamma_pi[1]*(1-gamma_pi[2]))
  
  
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
  
  tally <- tally + 1
  
  results_deterministic <- dplyr::bind_rows(results_deterministic,data.frame(alpha=a,beta=b,eps=E,T3=T,L=L2,M=M2))
}

# we may have more than 1000 rows for some reason
# but there are 1000 unique elements in alpha column
# we get the unique entries of our dataframe

results_deterministic <- results_deterministic %>%
  dplyr::filter(abs(alpha - lag(alpha)) > 1e-6 | is.na(lag(alpha))) %>%
  dplyr::ungroup()

write.csv(results_deterministic, "~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T3_alphabeta_deterministic.csv")
