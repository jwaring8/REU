rm(list=ls())
library(ggplot2)
library(pomp)

setwd("~/R/SimulatedData/SimulatedData")

## the data
# sample data was generated using this model and parameters below

S.0 <- 1e4-1
I.0 <- 1
N <- S.0+I.0
true_parameters <- c(beta_sd= 0.1, beta1 = 1.5, mu=0.1, gamma=0.6, rho= 0.8, 
                     N=N, S.0= S.0, I.0 = I.0, RealCases.0 = 1)

Fake<-read.csv(file='sample_data.csv', header=TRUE, sep=',')
data<-subset(Fake,select=c("time", "Cases"))
rm(Fake)

allSims <- list()

### I would have only time and Cases

### i added  stochasticity in births and white noise on the force of infection

## States: S I R, N observations are Cases, Time is time.

#beta sd 
#beta1
#mu
#gamma
#rho (observation probability)

## Writing the C Snippet


step.fun<- Csnippet("
                    double dQ = rgammawn(beta_sd,dt);
                    double lambda = beta1*I/N*dQ/dt;
                    double births = rpois(mu*N);
                    
                    S += (births - lambda*S - mu*S)*dt;
                    I += (lambda*S - gamma*I - mu*I)*dt;
                    RealCases += gamma*I*dt;
                    
                    ")
### create pomp object with the states and the paraeters data and process, obsnames is the observed cases and zeronames is Real cases


# Writing the rmeasure

rmeas <- Csnippet("
                  Cases = rpois(rho*RealCases + 1e-2);
                  ")


trans1 <- Csnippet("
                   Tbeta_sd = log(beta_sd);
                   Tbeta1 = log(beta1);
                   Tmu = log(mu);
                   Tgamma = log(gamma);
                   Trho = logit(rho);
                   TN = log(N);
                   TS_0 = log(S_0);
                   TI_0 = log(I_0);
                   TRealCases_0 = log(RealCases_0);
                   
                   ")

trans2 <- Csnippet("
                   Tbeta_sd = exp(beta_sd);
                   Tbeta1 = exp(beta1);
                   Tmu = exp(mu);
                   Tgamma = exp(gamma);
                   Trho = 1/(1+exp(-rho));
                   TN = exp(N);
                   TS_0 = exp(S_0);
                   TI_0 = exp(I_0);
                   TRealCases_0 = exp(RealCases_0);
                   
                   ")



po<-pomp(data=data, time="time", t0=0, 
         rprocess=euler.sim(step.fun=step.fun,delta.t=0.1), statenames= c("S", "I", "RealCases"),
         toEstimationScale=trans1,
         fromEstimationScale=trans2,
         rmeasure = rmeas,
         paramnames=c("beta_sd", "beta1", "mu", "gamma","rho", "N", "S.0", "I.0", "RealCases.0"), 
         obsnames= c("Cases"), 
         zeronames=c("RealCases"))


# Defining d measure we add a quantity to make sure it doesnt go to zero 

dmeas <- Csnippet("
                  lik = dpois(Cases,rho*RealCases+ 1e-2,give_log);
                  ")

po<-pomp(data=data, time="time", t0=0, 
         rprocess=euler.sim(step.fun=step.fun,delta.t=0.1), 
         statenames= c("S", "I", "RealCases"),
         rmeasure = rmeas,
         dmeasure=dmeas,
         toEstimationScale=trans1,
         fromEstimationScale=trans2,
         paramnames=c("beta_sd", "beta1", "mu", "gamma","rho", "N", "S.0", "I.0", "RealCases.0"), 
         obsnames= c("Cases"), 
         zeronames=c("RealCases"))

### Attempting a deterministic model using the skeleton

skel<-Csnippet("
double lambda = beta1*I/N;
double births = mu*N;

DS = births - lambda*S -mu*S;
DI = lambda*S -gamma*I -mu*I;
DRealCases = gamma*I;

")


po<-pomp(data=data, time="time", t0=0, 
         rprocess=euler.sim(step.fun=step.fun,delta.t=0.1), 
         skeleton=skel,skeleton.type="vectorfield",
         statenames= c("S", "I", "RealCases"),
         rmeasure = rmeas,
         dmeasure=dmeas,
         toEstimationScale=trans1,
         fromEstimationScale=trans2,
         paramnames=c("beta_sd", "beta1", "mu", "gamma","rho", "N", "S.0", "I.0", "RealCases.0"), 
         obsnames= c("Cases"), 
         zeronames=c("RealCases"))

for(i in 2:25){
  newSimulation <- simulate(po,nsim=1,params=c(beta_sd= 0.1, beta1 = 1.5, mu=0.1, gamma=0.6, rho= 0.8, N=N, S.0= S.0, I.0 = I.0, RealCases.0 = 1),
                            as.data.frame=T)
  
  ggplot(data=newSimulation, mapping=aes(x=time, y = Cases)) +
    geom_line()+geom_point()+expand_limits(y=0)+theme_classic()
  allSims[[length(allSims)+1]] <- newSimulation
}

save(allSims, file="simulations.Rda")


