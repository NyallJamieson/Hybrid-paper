# Second plot
# Ctmc sample of incubation times
MODEL.onestep <- function(x,params){
  L <- x[2]
  M <- x[3]
  alpha <- params["alpha"]
  beta <- params["beta"]
  lambda <- params["lambda"]
  rates <- c(
    sphago <- alpha*L,
    fphago <- beta*L,
    burst <- lambda*M
  )
  total.rates <- sum(rates)
  if (total.rates==0){
    tau <- -Inf
  }else{
    tau <- rexp(n=1,rate=total.rates)}
  transitions <- list(
    successfulPHA <- c(-1,1),
    failedPHA <- c(-1,0),
    bursted <- c(62,-1)
  )
  event <- sample.int(n=3,size=1,prob=rates/total.rates)
  x+c(tau,transitions[[event]])
}
illcount <- 0
illtime <- c()
varName1 <- "illcount"
varName2 <- "illtime"
MODEL.simul <- function(y,maxstep=10000000,threshold=50661){
  dose <- xstart["L"]
  x <- xstart
  illcount.copy <- get(varName1)
  illtime.copy <- get(varName2)
  output <- array(dim=c(maxstep+1,3))
  colnames(output) <- names(x)
  output[1,] <- x
  k <- 1
  while ((k <= maxstep) && (x["L"] > 0 | x["M"] > 0) && (x["L"] < threshold)){
    k <- k+1
    output[k,] <- x <- MODEL.onestep(x,params)}
  t <- output[k,1]
  if (x["L"]!=0){
    return(c(x["time"]))
  } else {return(c(NA))}
}
params <- c(alpha=0.08944449,beta=1.088102,lambda=0.04050593)
results <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(results) <- c("L","time")
xstart <- c(time=0,L=1,M=0)
results <- as.data.frame(na.omit(do.call(rbind,parallel::mclapply(1:2000,MODEL.simul,mc.cores=5)),na.action="omit"))
times_ctmc <- results[,1]

# Mean growth curve
thres <- 50661
a <- 0.089
b <- 1.088
h <- 0.041
g <- a+b
G <- 62.792
r <- 32.072
dt <- 0.01
model <- function(t,y,parms){
  dy1 <- -(a+b)*y[1]+h*G*y[2]
  dy2 <- a*y[1]-h*y[2]
  list(c(dy1,dy2))
}
yini <- c(y1=1,y2=0)
times <- seq(from=0,to=110,by=0.01)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)
times_ode <- out1[length(which(out1[,2]<thres))+1,1]

t <- seq(0,110,by=0.01)

# Sde simulation
simulate_sde <- function(y){
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
  }
  
  return(as.numeric(l)) # Return the full trajectories of l
}

times_sde <- c()
while (length(times_sde)<2000){
  print(length(times_sde))
  A <- simulate_sde(1)
  if (length(which(A>50661))>0){
    times_sde <- append(times_sde, t[length(which(A<50661))+1])
  }
}

ker_ctmc <- density(times_ctmc)
ker_sde <- density(times_sde,"bw"=1.5)

# Save plot with 600 DPI resolution
tiff("Fig2.2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(
  ker_ctmc,
  xlab = "Hours post infection",
  main = "Incubation-period distribution",
  lwd = 3,
  ylim = c(0, 0.1),
  col = "skyblue",
  cex.lab = 1.5,   # Larger font size for axis labels
  cex.main = 1.5,  # Larger font size for the main title
  cex.axis = 1.2   # Larger font size for axis tick labels
)
lines(ker_sde, lwd = 3, col = "red")
lines(c(times_ode, times_ode), c(0, 1), col = "green", lwd = 3)

# Increase legend text size using cex
legend(
  "topright",
  col = c("skyblue", "red", "green"),
  legend = c("MTBP", "SDE", "ODE"),
  lwd = 3,
  cex = 1.5 # Larger font size for the legend
)

# Close the device
dev.off()

# Third plot now
# Save plot with 600 DPI resolution
tiff("Fig2.3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(
  log = "y",
  out1[1:which(out1[, 1] == times_ode), 1],
  out1[1:which(out1[, 1] == times_ode), 2],
  lwd = 3,
  type = "l",
  xlab = "Hours post infection",
  ylab = "Extracellular Legionella population",
  main = "Extracellular Legionella population over time",
  xlim = c(0, 100),
  ylim = c(1e-02, 1e+06),
  cex.lab = 1.5,    # Larger font size for axis labels
  cex.main = 1.5,   # Larger font size for the main title
  cex.axis = 1.2    # Larger font size for axis ticks
)

for (i in 1:4) {
  A <- simulate_sde(1)
  lines(
    log = "y",
    out1[1:which(out1[, 1] == times_ode), 1],
    A[1:which(out1[, 1] == times_ode)],
    lwd = 1,
    col = i + 1
  )
}

lines(
  log = "y",
  out1[1:which(out1[, 1] == times_ode), 1],
  out1[1:which(out1[, 1] == times_ode), 2],
  lwd = 3,
  type = "l",
  xlab = "Hours post infection",
  ylab = "Extracellular Legionella population",
  main = "Extracellular Legionella population over time",
  xlim = c(0, 100),
  ylim = c(1e-02, 1e+06)
)

# Increase legend text size using cex
legend(
  "bottomright",
  col = c("black", 2, 3, 4, 5),
  legend = c("ODE model", "SDE trajectory 1", "SDE trajectory 2", "SDE trajectory 3", "SDE trajectory 4"),
  lwd = 3,
  cex = 1.5 # Larger font size for the legend
)

# Close the device
dev.off()

# Fourth plot now
# mean growth curve
thres <- 50661
a <- 0.089
b <- 1.088
h <- 0.041
g <- a+b
G <- 62.792
r <- 32.072
dt <- 0.01
model <- function(t,y,parms){
  dy1 <- -(a+b)*y[1]+h*G*y[2]
  dy2 <- a*y[1]-h*y[2]
  list(c(dy1,dy2))
}
yini <- c(y1=400,y2=0)
times <- seq(from=0,to=110,by=0.01)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)

t <- seq(0,110,by=0.01)
# Sde simulation
simulate_sde <- function(y){
  l <- rep(0, length(t))
  m <- rep(0, length(t))
  l[1] <- 400
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
  }
  
  return(as.numeric(l)) # Return the full trajectories of l
}

# Save plot with 600 DPI resolution
tiff("Fig2.4.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(
  log = "y",
  out1[1:which(out1[, 1] == times_ode), 1],
  out1[1:which(out1[, 1] == times_ode), 2],
  lwd = 3,
  type = "l",
  xlab = "Hours post infection",
  ylab = "Extracellular Legionella population",
  main = "Extracellular Legionella population over time",
  xlim = c(0, 100),
  ylim = c(1e+0, 1e+8),
  cex.lab = 1.5,    # Larger font size for axis labels
  cex.main = 1.5,   # Larger font size for the main title
  cex.axis = 1.2    # Larger font size for axis ticks
)

for (i in 1:4) {
  A <- simulate_sde(1)
  lines(
    log = "y",
    out1[1:which(out1[, 1] == times_ode), 1],
    A[1:which(out1[, 1] == times_ode)],
    lwd = 1,
    col = i + 1
  )
}

lines(
  log = "y",
  out1[1:which(out1[, 1] == times_ode), 1],
  out1[1:which(out1[, 1] == times_ode), 2],
  lwd = 3,
  type = "l"
)

# Increase legend text size using cex
legend(
  "bottomright",
  col = c("black", 2, 3, 4, 5),
  legend = c("ODE model", "SDE trajectory 1", "SDE trajectory 2", "SDE trajectory 3", "SDE trajectory 4"),
  lwd = 3,
  cex = 1.5 # Larger font size for the legend
)

# Close the device
dev.off()

# First plot now
MODEL.simul2 <- function(y,maxstep=10000000,threshold=50661){
  done <- 0
  while (done==0){
    dose <- xstart["L"]
    x <- xstart
    illcount.copy <- get(varName1)
    illtime.copy <- get(varName2)
    output <- array(dim=c(maxstep+1,3))
    colnames(output) <- names(x)
    output[1,] <- x
    k <- 1
    while ((k <= maxstep) && (x["L"] > 0 | x["M"] > 0) && (x["L"] < threshold)){
      k <- k+1
      output[k,] <- x <- MODEL.onestep(x,params)}
    t <- output[k,1]
    if (x["L"]!=0){
      done <- 1
    }
  }
  return(output)
}

simulate_sde_3 <- function(y){
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
  }
  
  return(as.numeric(l)) # Return the full trajectories of l
}

#A2 <- MODEL.simul2(1)
#saveRDS(A2,"MotivationBP.rds")
A2 <- readRDS("~/PhD-work/Chapter 3/Data/Figure 2A motivation/MotivationBP.rds")
#B <- simulate_sde_3(1)
#write.csv(B,"MotivationSDE.csv")
B <- read.csv("~/PhD-work/Chapter 3/Data/Figure 2A motivation/MotivationSDE.csv")[,2]

# Save plot with 600 DPI resolution
tiff("Fig2.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(
  A2[which(A2[, 1] < 25), 1],
  A2[which(A2[, 1] < 25), 2],
  type = "l",
  xlab = "Hours post infection",
  ylab = "Extracellular Legionella population",
  main = "Extracellular Legionella population over time",
  lwd = 3,
  cex.lab = 1.5,   # Larger font size for axis labels
  cex.main = 1.5,  # Larger font size for the main title
  cex.axis = 1.2   # Larger font size for axis ticks
)

lines(out1[1:2501, 1], B[1:2501], col = "red", lwd = 2)
lines(A2[which(A2[, 1] < 25), 1], A2[which(A2[, 1] < 25), 2], col = "black", lwd = 2)

# Increase legend text size using cex
legend(
  "topleft",
  col = c("black", "red"),
  legend = c("Markov chain trajectory", "SDE trajectory"),
  lwd = 3,
  cex = 1.5 # Larger font size for the legend
)

# Close the device
dev.off()
