#This file will generate a file called "allRawTrajMatches.Rda" containing trajectory matches 
#from 25 simulations of simple epidemic data with the following parameters: beta_sd = 0.1,
#Beta=1.5, mu=0.1, gamma=0.6, rho=0.8, N=10000. This also uses raw incidence data.

require(pomp)
require(ggplot2)
require(magrittr)

setwd("~/R/SimulatedData/SimulatedData")

#Loading Rda file with all the simulations
load("simulations.Rda") #generated from Simulator.R

allRawTrajMatches <- list() #List to keep track of the 39 traj matches for each of the 25 simulations 

for(i in 1:25){
  if(i ==1) {
    #Read in sample data file and plot the epidemic curve
    sampleData <- read.csv(file ="sample_data.csv")
    ggplot(data=sampleData, mapping=aes(x=time, y= Cases)) +
      geom_line()+geom_point()+expand_limits(y=0)+theme_classic()
    sampleData <- subset(sampleData, select=c("time", "Cases"))
  } else {
    #Read in sample data file and plot the epidemic curve
    sampleData <- allSims[[i-1]]
    ggplot(data=sampleData, mapping=aes(x=time, y= Cases)) +
      geom_line()+geom_point()+expand_limits(y=0)+theme_classic()
    sampleData <- subset(sampleData, select=c("time", "Cases"))
   }
  #Continuous-time process model
  step.fun <- Csnippet("
                       double dQ = rgammawn(beta_sd,dt);
                       double lambda = Beta*I/N*dQ/dt;
                       double births = rpois(mu*N);
                       
                       S += (births - lambda*S - mu*S)*dt;
                       I += (lambda*S - gamma*I - mu*I)*dt;
                       RealCases += gamma*I*dt;          
                       ")
  
  #Initializer
  init <- Csnippet("
                   S = N-1;
                   I= 1;
                   RealCases = 1;
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
  #CIs <- list()
  
  for(i in 2:nrow(sampleData)) {
    #Construction of pomp model 
    sample <- pomp(sampleData[1:i,], time="time", t0=0, 
                   rprocess=euler.sim(step.fun, delta.t=0.1), initializer = init,
                   statenames=c("S", "I", "RealCases"), paramnames=c("Beta", "gamma", "mu", "rho", "N","beta_sd"),
                   obsnames = "Cases", zeronames = "RealCases", rmeasure = rmeas, dmeasure = dmeas, 
                   skeleton = vectorfield(skel), toEstimationScale=logtrans, fromEstimationScale = exptrans)
    
    #Trajectory Matching 
    tm <- traj.match(sample, start=c(beta_sd=0.12, Beta=1.8, mu=0.15, gamma=0.65, rho=0.75, N=10000), 
                     est=c("beta_sd", "Beta", "mu", "gamma", "rho", "N"), transform = TRUE)
    #Simulate fitted model against data
    p <- simulate(tm,nsim=10,as.data.frame=TRUE,include.data=TRUE)
    p <- subset(p, select=c("time", "Cases", "sim"))
    p <- subset(p, !is.na(p$Cases))
    ggplot(data=p,aes(x=time,y=Cases,group=sim,alpha=(sim=="data")))+
      scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                         labels=c(`FALSE`="simulation",`TRUE`="data"))+
      geom_line() -> newPlot
    trajMatches[[length(trajMatches)+1]] <- newPlot
    #CIs[[length(CIs)+1]] <- confint(tm, parm = c("beta_sd", "Beta", "mu", "gamma", "rho", "N"), level = 0.95)
    
  } #for loop
  
  allRawTrajMatches[[length(allRawTrajMatches)+1]] <- trajMatches
  
  #Testing model with simulations
  sims <- simulate(sample,params=c(beta_sd=0.1, Beta=1.5, mu=0.1, gamma=0.6, rho=0.8, N=10000),
                   nsim=20,as.data.frame=TRUE,include.data=TRUE)
  
  ggplot(sims,mapping=aes(x=time,y=Cases,group=sim,color=sim=="data"))+
    geom_line()+guides(color=FALSE)
  
}

save(allRawTrajMatches, file="allRawTrajMatches.Rda")