setwd("~/Documents/Jonathan's Cholera Work")

library(XLConnect)
library(ggplot2)
library(pomp)
library(magrittr)

#Read in data from file and plot epidemic curve
regionalHaitiCholeraData <- readWorksheet(loadWorkbook("Haiti_cholera-cases-per-month-and-per-region-since-2015.xlsx"), sheet = 1)
choleraHaitiList <- split(regionalHaitiCholeraData, regionalHaitiCholeraData$Departement)
ggplot(data=choleraHaitiList[[1]], mapping=aes(x=choleraHaitiList[[1]]$Month, y=choleraHaitiList[[1]]$Cases)) + 
  geom_line() + geom_point() + expand_limits(y=0) + theme_classic() + xlab("Month") + 
  ylab("Cases") + ggtitle("Monthly Cholera Cases in the Artibonite Department of Haiti: January 2015-March 2016")

#Continuous Time Process Model
sir_step <- Csnippet(" 
  double lambda, beta;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection

  // Transitions
  // From class S
  double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections
  // From class I
  double transI = rbinom(I, 1.0 - exp(- gamma * dt)); // No of transitions I->R
  // Balance the equations
  S -= transS;
  I += transS - transI;
  R += transI;
  H += transI; //track incidence 
")

#Initializer
init <- Csnippet("
  S = N - 391;
  I = 391;
  R = 0;
  H = 0;
")

#Adjusting time for pomp model (Converting POSIXct to numeric time)
choleraHaitiList[[1]]$Month <- replace(as.numeric(choleraHaitiList[[1]]$Month), values = 1:15)
choleraHaitiList[[1]]$Month

#Adjusting data to be used
sampleData <- subset(choleraHaitiList[[1]], select = c("Cases", "Month"))

#First construction of pomp model
pompModel <- pomp(data=sampleData, time = "Month", t0=0, rprocess=euler.sim(sir_step, delta.t=1/30), 
                  initializer = init, paramnames = c("N", "R0", "gamma"), statenames = c("S", "I", "R", "H"), obsnames = "Cases")

#Rmeasure model
rmeas <- Csnippet("
  Cases = rpois(rho * H);                  
")

#Adding rmeasure to pomp model
pompModel <- pomp(pompModel, rmeasure=rmeas, statenames="H", paramnames ="rho", obsnames = "Cases")

#Dmeasure model 
dmeas <- Csnippet("
  lik = dpois(Cases, H * rho, give_log);
")

#Adding dmeasure to pomp model
pompModel <- pomp(pompModel, dmeasure=dmeas, statenames="H", paramnames="rho", obsnames="Cases")

#Testing model with simulations
sims <- simulate(pompModel,params=c(gamma=0.1, rho=0.8, N=1000, R0 = 1.5),
                 nsim=20,as.data.frame=TRUE,include.data=TRUE)

ggplot(sims,mapping=aes(x=time,y=Cases,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

#Deterministic Model Skeleton
skel <- Csnippet("
  double lambda, beta;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection

  // Balance the equations
  DS = - lambda * S;
  DI = (lambda * S) - (gamma * I);
  DR = gamma * I;
  DH = gamma * I;
")

#Adding skeleton to pomp model
pompModel <- pomp(pompModel, skeleton=vectorfield(skel), statenames=c("S", "I", "R", "H"), 
               paramnames=c("gamma","rho", "N", "R0"), obsnames="Cases")

#Paramter transformations
logtrans <- Csnippet("
Tgamma = log(gamma);
Trho = log(rho);
TN = log(N);
TR0 = log(R0);
")

exptrans <- Csnippet("
Tgamma = exp(gamma);
Trho = exp(rho);
TN = exp(N);
TR0 = exp(R0);
")

#Adding parameter transformations to pomp model
pompModel <- pomp(pompModel, toEstimationScale=logtrans, fromEstimationScale=exptrans,
               paramnames=c("gamma", "rho", "N", "R0"))

#Particle Filtering
pf <- pfilter(pompModel, Np=100000,params=c(gamma=0.6, rho=0.8, N=1000000, R0 = 1.5), filter.traj=TRUE)
logLik(pf)
coef(pf)


#Trajectory Matching 
tm <- traj.match(pf, est=c("gamma", "rho", "N", "R0"), transform = TRUE)
summary(tm)

#Simulate fitted model against data
simulate(tm,nsim=10,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=Cases,group=sim,alpha=(sim=="data")))+
  scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                     labels=c(`FALSE`="simulation",`TRUE`="data"))+
  geom_line()

