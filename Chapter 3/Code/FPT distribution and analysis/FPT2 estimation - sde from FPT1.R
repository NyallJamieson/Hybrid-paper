# Set parameter values
t <- seq(0,110,by=0.01)
dt <- 0.01
a <- 0.089
b <- 1.088
h <- 0.041
m <- 62.792
G <- m
r <- 32.072

# Sde simulation
simulate_sde <- function(y){
  
  done <- NULL
  
  l <- rep(0, length(t))
  m <- rep(0, length(t))
  l[1] <- 49.12922
  m[1] <- 25.01429
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
    
    if (l[i]>92.1405 && m[i]>46.91361){
      time <- t[i]
      done <- time
      break
    }
  }
  return(done) # Return the time from FPT1 to FPT2
}

# Run simulation 10000 times
num_sim <- 10000
num_cores <- 2
results <- as.numeric(parallel::mclapply(1:num_sim, simulate_sde, mc.cores = num_cores))

# Read in T1 data from Gillespie algorithm
T1_gill <- read.csv("~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T1_gill.csv")[,2]

# Take convolution of this sde with the Gillespie data
FPT2_sde_from_FPT1 <- c()
for (i in 1:10000){
  print(i)
  FPT2_sde_from_FPT1 <- append(FPT2_sde_from_FPT1,sample(T1_gill,size=1,replace=TRUE)+sample(results,size=1,replace=TRUE))
}

# Now save
write.csv(FPT2_sde_from_FPT1,"~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T2_sde.csv")