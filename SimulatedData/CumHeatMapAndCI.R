#Generate heat map and confidence intervals for raw data simulations. 

require(pomp)
require(ggplot2)
require(plotly)
require(RColorBrewer)


sampleData <- read.csv(file ="sample_data.csv")
gamma <- 0.6
cumulativeSampleData <- sampleData
actualCumulative <- trunc(cumsum(cumulativeSampleData$Cases-(gamma*cumulativeSampleData$Cases)))
cumulativeSampleData <- cbind(cumulativeSampleData,actualCumulative)
cumulativeSampleData <- subset(cumulativeSampleData, select=c("time", "Cases", "actualCumulative"))
ggplot(data=cumulativeSampleData, mapping=aes(x=time, y= actualCumulative)) +
  geom_line()+geom_point()+expand_limits(y=0)+theme_classic()

#Continuous-time process model
step.fun <- Csnippet("
                     double dQ = rgammawn(beta_sd,dt);
                     double lambda = Beta*I/N*dQ/dt;
                     double births = rpois(mu*N);
                     
                     S += (births - lambda*S - mu*S)*dt;
                     I += (lambda*S - gamma*I - mu*I)*dt;
                     RealCases += gamma*I*dt;
                     CumulativeCases += (RealCases - (RealCases * gamma))*dt;
                     ")

#Initializer
init <- Csnippet("
                 S = N-1;
                 I= 1;
                 RealCases = 1;
                 CumulativeCases = 1;
                 ")

#Rmeasure model
rmeas <- Csnippet("
                  Cases = rpois(RealCases * rho);                
                  ")

#Dmeasure model
dmeas <- Csnippet("
                  lik = dpois(Cases, RealCases * rho, give_log);                  
                  ")
  
  #Deterministic model skeleton
  skel <- Csnippet("
                   double lambda = Beta*I/N;
                   double births = mu*N;
                   
                   DS = births - lambda*S -mu*S;
                   DI = lambda*S -gamma*I -mu*I;
                   DRealCases = gamma*I;    
                   DCumulativeCases = (RealCases - (RealCases * gamma));
                   ")

  #Paramter transformations
  logtrans <- Csnippet("
                       Tbeta_sd = log(beta_sd);
                       TBeta = log(Beta);
                       Tmu = log(mu);
                       Tgamma = log(gamma);
                       Trho = log(rho);
                       TN = log(N);                     
                       ")

  exptrans <- Csnippet("
                       Tbeta_sd = exp(beta_sd);
                       TBeta = exp(Beta);
                       Tmu = exp(mu);
                       Tgamma = exp(gamma);
                       Trho = exp(rho);
                       TN = exp(N);                     
                       ")
  
  sample <- pomp(cumulativeSampleData, time="time", t0=0, 
                 rprocess=euler.sim(step.fun, delta.t=0.1), initializer = init,
                 statenames=c("S", "I", "RealCases", "CumulativeCases"), paramnames=c("Beta", "gamma", "mu", "rho", "N","beta_sd"),
                 obsnames = "Cases", zeronames = "RealCases", rmeasure = rmeas, dmeasure = dmeas, 
                 skeleton = vectorfield(skel), toEstimationScale=logtrans, fromEstimationScale = exptrans)

#Vectors and matrix for looping over different beta and gamma values for log lik estimates
betaVector <- seq(from=1.0, to=2.0, by=((2.0-1.0)/99))
gammaVector <- seq(from=0.25, to=1.0, by=((1.0-.25)/99))
logLiks <- matrix(nrow = 100, ncol = 100)
  
#Filling logLiks matrix at different Beta and Gamma values 
if(file.exists("CumLogLiks.Rda")) {    
  load("CumLogLiks.Rda")
} else {
  for(i in 1:100) {
    for(j in 1:100) {
      sims <- simulate(sample,params=c(beta_sd=0.1, Beta=betaVector[i], mu=0.1, gamma=gammaVector[j], rho=0.8, N=10000),
                       nsim=500, states=TRUE)
      ll <- dmeasure(sample, y=obs(sample), x=sims, times=time(sample), log=TRUE, 
                     params=c(beta_sd=0.1, Beta=betaVector[i], mu=0.1, gamma=gammaVector[j], rho=0.8, N=10000))
      ll <- apply(ll, 1, sum)
      ll <- max(ll)
      logLiks[i,j] <- ll
    }
  }
    save(logLiks, file="CumLogLiks.Rda")
}

#Plotting logLiks matrix as a HeatMap
my_palette <- colorRampPalette(c("red", "yellow", "green", "blue", "purple"))(n = 1999)
logLikliehood <- matrix(nrow=100, ncol=100)
for(i in 1:100) {
  for(j in 1:100) {
    if(logLiks[i,j] >= -1000) {
      logLikliehood[i,j] = logLiks[i,j]
    }else{
      logLikliehood[i,j] = -1000
    }
  }
}
x <- list(
  title = "Beta",
  autotick = FALSE,
  ticks = "outside",
  tick0 = 0,
  dtick = 0.25
)
y <- list(
  title = "Gamma",
  autotick = FALSE,
  ticks = "outside",
  tick0 = 0,
  dtick = 0.25
)
plot_ly(x=betaVector, y=gammaVector, z=logLikliehood, type="contour", colors=my_palette) %>%
  layout(xaxis=x, yaxis=y)

#Plotting logLiks as a function of Beta and Gamma
vec <- vector(length=100)
for(i in 1:100){
  vec[i] <- logLiks[i,i]
}
y <- list(
  title = "Log Likliehood",
  autotick = FALSE,
  ticks = "outside",
  tick0 = 0,
  dtick = 250
)

plot_ly(x=betaVector, y=vec) %>%
  layout(xaxis=x, yaxis=y)
x <- list(
  title = "Gamma",
  autotick = FALSE,
  ticks = "outside",
  tick0 = 0,
  dtick = 0.25
)
y <- list(
  title = "Log Likelihood",
  autotick = FALSE,
  ticks = "outside",
  tick0 = 0,
  dtick = 250
)
plot_ly(x=gammaVector, y=vec) %>%
  layout(xaxis=x, yaxis=y)