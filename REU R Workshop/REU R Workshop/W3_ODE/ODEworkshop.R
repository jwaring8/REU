polynomial <- function(x) {
  y <- (5*(x^3))+(6*x)-(4/x)
  return(y)
}
polynomial(3)
x <- seq(-10, -1, .01)
y <- polynomial(x)

plot(x, y, type ='l', xlab = 'X values', ylab = 'Y values', main = '5x^3+6x-4/x', lwd=3.0)
x[which(y==-7)]

discreteTimeModel <- function(lambda, popAtTime) {
  return(lambda * popAtTime)
}

discreteTimeModel2 <- function(lambda, initialPop, t) {
  return(lambda^t * initialPop)
}

x2 <- 0:20
y2 <- discreteTimeModel2(0.75, 1, x2)
y3 <- discreteTimeModel2(1, 1, x2)
y4 <- discreteTimeModel2(2.0, 1, x2)
plot(x2, y2, pch = 21, bg = "grey40", col = 'black', cex = 2.0, xlab = 'X values', ylab ='Y values')
plot(x2, y3, col = 'black', pch = 21, bg = "dodgerblue", cex = 2.0)
plot(x2, y4, col = 'black', pch = 21, bg = "firebrick", cex = 2.0)

continuousTimeModel <- function(r, t, initialPop) {
  return(initialPop * 2.71828^(r*t))
}

y5 <- continuousTimeModel(0.05, x2, 1)

plot(x2, y5, col = 'black', pch =21, bg = "green", cex =2.0)

#LSODA stuff
y0 <- c(N = 1.0) #vector of initial variable values in the system (initialPop = 1)
time.vec <- seq(0, 10, 0.1) #vector of times for function to return the values of variables in the Differential equations
r = .75
params <- c(r = r) #vector of param values
func <- function(t, y, params) { #function to represent right hand side of diff equation
  N <- y["N"] #N is the value of "N" in the y vector
  r <- params["r"] #r is the value of "r" in the params vector
  dN <- r*N #dN/dt
  return(list(c(dN)))
}
odeoutput <- lsoda(y0, time.vec, func, params)
odeoutput <- as.data.frame(odeoutput)
plot(N~time, data=odeoutput, pch = 21, bg ='dodgerblue', col = 'black', cex=2.0)

#SIR model stuff
ClosedSIR <- function(t, y, params) {
  S <- y["S"]; I <- y["I"]; R<-y["R"]; N <- S+I+R
  beta <- params["beta"]; gamma <- params["gamma"]
  dS <- -beta*I/N*S
  dI <- beta*I/N*S - gamma*I
  dR <- gamma*I
  return(list(c(dS, dI, dR)))
}
t.vec <- seq(0, 120, 0.1)
beta <- 0.3; gamma <- 1/7
params <- c(beta = beta, gamma = gamma)
N <- 10000
I0 <- 1/N
S0 <- 9999/N
R0 <- 0/N
y0 <- c(S=S0, I=I0, R=R0)

SIRoutput <- lsoda(y0, t.vec, ClosedSIR, params)
SIRoutput <- as.data.frame(SIRoutput)
plot(S~time, data=SIRoutput, type='b', col='red', ylim=c(0,1), ylab="Proportion", xlab='Time')
lines(I~time, data=SIRoutput, type='b', col = 'blue')
lines(R~time, data=SIRoutput, type='b', col='green')

#Reed-Frost chain binomial model stuff
alpha <- 0.96
N <- 100
i <- 1
s0 <- N-i 
t <- 0 #start time
out.rf <- c(t, s0, i)
while(i > 0){
  s <- rbinom(1, s0, alpha^i)
  i <- s0-s
  t <- t+1
  out.rf <- rbind(out.rf, c(t,s,i))
  s0 <- s
}

out.rf <- as.data.frame(out.rf)
names(out.rf) <- c("time", "S", "I")

#SIS Model
ClosedSIS <- function(t, y, params) {
  S <- y["S"]; I <- y["I"]; N <- S+I
  beta <- params["beta"]; gamma <- params["gamma"]
  dS <- -beta*I/N*S + gamma*I
  dI <- beta*I/N*S - gamma*I
  return(list(c(dS, dI)))
}
t.vec <- seq(0, 120, 0.1)
beta <- 0.4; gamma <- 1/7
params <- c(beta = beta, gamma = gamma)
N <- 10000
I0 <- 1/N
S0 <- 9999/N
y0 <- c(S=S0, I=I0)

SISoutput <- lsoda(y0, t.vec, ClosedSIS, params)
SISoutput <- as.data.frame(SISoutput)
plot(S~time, data=SISoutput, type='b', col='red', ylim=c(0,1), ylab="Proportion", xlab='Time')
lines(I~time, data=SISoutput, type='b', col = 'blue')
