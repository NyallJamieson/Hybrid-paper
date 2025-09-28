# Set initial parameters

T2 <- 45.1
a <- 0.089
b <- 1.088
h <- 0.041
m <- 62.792
G <- m
r <- 32.072

# Mean
model <- function(t,y,parms){
dy1 <- -(a+b)*y[1]+h*G*y[2]
dy2 <- a*y[1]-h*y[2]
list(c(dy1,dy2))
}
yini <- c(y1=1,y2=0)
times <- seq(from=0,to=150,by=0.1)
out1 <- deSolve::ode(times=times,y=yini,func=model,parms=NULL)

L2 <- out1[length(which(times<=T2)),2]
M2 <- out1[length(which(times<=T2)),3]

# Tau leap simulation for Jump process
times2 <- c()

while (length(times2)<10000){

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
      times2 <- append(times2,t[i])
      print(length(times2))
      break
      }
    }
  }

# Now save
write.csv(times2,"~/PhD-work/Chapter 3/Data/FPT dataset for distribution/T1_tau.csv")
