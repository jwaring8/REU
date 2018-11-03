#This file will generate a file called "allCumTrajMatches.Rda" containing trajectory matches 
#from 25 simulations of simple epidemic data with the following parameters: beta_sd = 0.1,
#Beta=1.5, mu=0.1, gamma=0.6, rho=0.8, N=10000. This also uses raw incidence data.

require(pomp)
require(ggplot2)
require(magrittr)
require(plyr)

setwd("~/R/SimulatedData/SimulatedData")

#Loading Rda file with all the simulations
load("simulations.Rda") #generated from Simulator.R

allCumTrajMatches <- list()
for(i in 1:25){
  if(i ==1) {
    #Read in sample data file and plot the epidemic curve
    sampleData <- read.csv(file ="sample_data.csv")
    gamma <- 0.6
    cumulativeSampleData <- sampleData
    actualCumulative <- trunc(cumsum(cumulativeSampleData$Cases-(gamma*cumulativeSampleData$Cases)))
    cumulativeSampleData <- cbind(cumulativeSampleData,actualCumulative)
    cumulativeSampleData <- subset(cumulativeSampleData, select=c("time", "Cases", "actualCumulative"))
    ggplot(data=cumulativeSampleData, mapping=aes(x=time, y= actualCumulative)) +
      geom_line()+geom_point()+expand_limits(y=0)+theme_classic()
  } else {
    #Read in sample data file and plot the epidemic curve
    sampleData <- allSims[[i-1]]
    gamma <- 0.6
    cumulativeSampleData <- sampleData
    actualCumulative <- trunc(cumsum(cumulativeSampleData$Cases-(gamma*cumulativeSampleData$Cases)))
    cumulativeSampleData <- cbind(cumulativeSampleData,actualCumulative)
    cumulativeSampleData <- subset(cumulativeSampleData, select=c("time", "Cases", "actualCumulative"))
    ggplot(data=cumulativeSampleData, mapping=aes(x=time, y= actualCumulative)) +
      geom_line()+geom_point()+expand_limits(y=0)+theme_classic()
  }
  
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
  
  trajMatches <- list()
  
  for(j in 2:nrow(cumulativeSampleData)) {
    #Construction of pomp model 
    sample <- pomp(cumulativeSampleData[1:j,], time="time", t0=0, 
                   rprocess=euler.sim(step.fun, delta.t=0.1), initializer = init,
                   statenames=c("S", "I", "RealCases", "CumulativeCases"), paramnames=c("Beta", "gamma", "mu", "rho", "N","beta_sd"),
                   obsnames = "Cases", zeronames = "RealCases", rmeasure = rmeas, dmeasure = dmeas, 
                   skeleton = vectorfield(skel), toEstimationScale=logtrans, fromEstimationScale = exptrans)
    #Trajectory Matching 
    tm <- traj.match(sample, start=c(beta_sd=0.12, Beta=1.8, mu=0.15, gamma=0.65, rho=0.75, N=10000), 
                     est=c("beta_sd", "Beta", "mu", "gamma", "rho", "N"), transform = TRUE)
    #Simulate fitted model against data
    p <- simulate(tm,nsim=10,as.data.frame=TRUE,include.data=TRUE)
    p$CumulativeCases[1:j] <- p$actualCumulative[1:j]
    p <- subset(p, select=c("time", "CumulativeCases", "sim"))
    p <- subset(p, !is.na(p$CumulativeCases))
    ggplot(data=p,aes(x=time,y=CumulativeCases,group=sim,alpha=(sim=="data")))+
      scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                         labels=c(`FALSE`="simulation",`TRUE`="data"))+
      geom_line() -> newPlot
     trajMatches[[length(trajMatches)+1]] <- newPlot
  }
  
  allCumTrajMatches[[length(allCumTrajMatches)+1]] <- trajMatches
  
}

save(allCumTrajMatches, file="allCumTrajMatches.Rda")