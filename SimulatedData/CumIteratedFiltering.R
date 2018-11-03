#Cumulative Iterated Filtering

require(pomp)
require(ggplot2)
require(magrittr)
require(foreach)
require(iterators)
require(doParallel)
require(plyr)
require(coda)
require(reshape2)

#Read in sample data file and plot the epidemic curve
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

#Dprior for posterior sampling
dprior <- Csnippet("
    lik = dunif(beta_sd, 0.01, 1.0, 1)+dunif(Beta, 1.0, 3.0, 1)+dunif(mu, 0.01, 1.0, 1)+
          dunif(gamma, 0.1, 2.0, 1)+dunif(rho, 0.05, 0.99, 1)+dunif(N, 1000, 10000, 1);
    lik = (give_log) ? lik : exp(lik);
    ")

sample <- pomp(cumulativeSampleData, time="time", t0=0, 
               rprocess=euler.sim(step.fun, delta.t=0.1), initializer = init,
               statenames=c("S", "I", "RealCases", "CumulativeCases"), paramnames=c("Beta", "gamma", "mu", "rho", "N","beta_sd"),
               obsnames = "Cases", zeronames = "RealCases", rmeasure = rmeas, dmeasure = dmeas, 
               skeleton = vectorfield(skel), dprior = dprior, 
               toEstimationScale=logtrans, fromEstimationScale = exptrans)

#Testing model with simulations
sims <- simulate(sample,params=c(beta_sd=0.12, Beta=1.8, mu=0.15, gamma=0.65, rho=0.75, N=10000),
                 nsim=20,as.data.frame=TRUE,include.data=TRUE)
sims$CumulativeCases[1:39] <- sims$actualCumulative[1:39]

ggplot(sims,mapping=aes(x=time,y=CumulativeCases,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

#IF2 algorithm
registerDoParallel()
bake(file="CumFakeData-mif.rds", {
  guesses <- sobolDesign(lower=c(beta_sd=0.05, Beta=1.0, mu=0.01, gamma=0.1, rho=0.05, N=10000),
                         upper=c(beta_sd=0.5, Beta=3.0, mu=1.0, gamma=0.99, rho=0.99, N=10000),
                         nseq=100)
  foreach (guess=iter(guesses,"row"),.combine=rbind,
           .options.mpi=list(seed=334065675),
           .packages=c("pomp","magrittr"),.errorhandling="remove") %dopar% {
             sample %>%
               mif2(start=unlist(guess),Nmif=50,Np=1000,transform=TRUE,
                    cooling.fraction.50=0.8,cooling.type="geometric",
                    rw.sd=rw.sd(beta_sd=0.02, Beta=0.02, mu=0.02, gamma=0.02, rho=0.02, N=0.02)) %>%
               mif2() -> mf
             ll <- logmeanexp(replicate(5,logLik(pfilter(mf))),se=TRUE)
             data.frame(loglik=ll[1],loglik.se=ll[2],as.list(coef(mf)))
           }
}) -> mles

mles %>% 
  subset(loglik==max(loglik)) %>% unlist() -> mle
simulate(sample,params=mle,nsim=10,as.data.frame=TRUE,include.data=TRUE) -> sims2
sims2$CumulativeCases[1:39] <- sims2$actualCumulative[1:39]
  ggplot(data=sims2, mapping=aes(x=time,y=CumulativeCases,group=sim,alpha=(sim=="data")))+
  scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                     labels=c(`FALSE`="simulation",`TRUE`="data"))+
  geom_line()

#Pairwise scatterplot matrix of MLES
pairs(~loglik+beta_sd+Beta+mu+gamma+rho+N,data=mles)

#Constructing profile likelihood for Beta
profileDesign(
  Beta=seq(from=1.0, to=3.0, length=20),
  lower=c(beta_sd=0.05, mu=0.01, gamma=0.1, rho=0.05, N=10000),
  upper=c(beta_sd=0.5, mu=1.0, gamma=0.99, rho=0.99, N=10000),
  nprof=50
) -> betaPD

bake("CumFakedata_Beta-profile.rds",{
  foreach (p=iter(betaPD,"row"),
           .combine=rbind,
           .errorhandling="remove",
           .packages=c("pomp","magrittr","reshape2","plyr"),
           .inorder=FALSE,
           .options.mpi=list(seed=1680158025)
  ) %dopar% {
    sample %>% 
      mif2(start=unlist(p),Nmif=50,Np=1000,transform=TRUE,
           cooling.fraction.50=0.8,cooling.type="geometric",
           rw.sd=rw.sd(beta_sd=0.02, mu=0.02, gamma=0.02, rho=0.02, N=0.02)) %>%
      mif2() -> mf
    
    pf <- replicate(5,pfilter(mf,Np=1000))  ## independent particle filters
    ll <- sapply(pf,logLik)
    ll <- logmeanexp(ll,se=TRUE)
    nfail <- sapply(pf,getElement,"nfail")  ## number of filtering failures
    
    data.frame(as.list(coef(mf)),
               loglik = ll[1],
               loglik.se = ll[2],
               nfail.min = min(nfail),
               nfail.max = max(nfail))
  } %>% arrange(Beta,-loglik)
}) -> beta_prof

pairs(~loglik+Beta+beta_sd+mu+gamma+rho+N,data=beta_prof,subset=loglik>max(loglik)-10)

beta_prof %>% 
  mutate(Beta=signif(Beta,8)) %>%
  ddply(~Beta,subset,loglik==max(loglik)) %>%
  ggplot(aes(x=Beta,y=loglik))+geom_point()+geom_smooth()

beta_prof %>% 
  subset(loglik==max(loglik)) %>% unlist() -> mle
simulate(sample,params=mle,nsim=10,as.data.frame=TRUE,include.data=TRUE) -> betaSim
betaSim$CumulativeCases[1:39] <- betaSim$actualCumulative[1:39]
ggplot(data=betaSim,mapping=aes(x=time,y=CumulativeCases,group=sim,alpha=(sim=="data")))+
scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                   labels=c(`FALSE`="simulation",`TRUE`="data"))+
geom_line()

#Constructing profile likelihood for gamma
profileDesign(
  gamma=seq(from=0.1, to=1.0, length=20),
  lower=c(beta_sd=0.05, mu=0.01, Beta=1.0, rho=0.05, N=10000),
  upper=c(beta_sd=0.5, mu=1.0, Beta=3.0, rho=0.99, N=10000),
  nprof=50
) -> gammaPD

bake("CumFakedata_gamma-profile.rds",{
  foreach (p=iter(gammaPD,"row"),
           .combine=rbind,
           .errorhandling="remove",
           .packages=c("pomp","magrittr","reshape2","plyr"),
           .inorder=FALSE,
           .options.mpi=list(seed=1680158025)
  ) %dopar% {
    sample %>% 
      mif2(start=unlist(p),Nmif=50,Np=1000,transform=TRUE,
           cooling.fraction.50=0.8,cooling.type="geometric",
           rw.sd=rw.sd(beta_sd=0.02, mu=0.02, Beta=0.02, rho=0.02, N=0.02)) %>%
      mif2() -> mf
    
    pf <- replicate(5,pfilter(mf,Np=1000))  ## independent particle filters
    ll <- sapply(pf,logLik)
    ll <- logmeanexp(ll,se=TRUE)
    nfail <- sapply(pf,getElement,"nfail")  ## number of filtering failures
    
    data.frame(as.list(coef(mf)),
               loglik = ll[1],
               loglik.se = ll[2],
               nfail.min = min(nfail),
               nfail.max = max(nfail))
  } %>% arrange(gamma,-loglik)
}) -> gamma_prof

pairs(~loglik+gamma+beta_sd+mu+Beta+rho+N,data=beta_prof,subset=loglik>max(loglik)-10)

gamma_prof %>% 
  mutate(gamma=signif(gamma,8)) %>%
  ddply(~gamma,subset,loglik==max(loglik)) %>%
  ggplot(aes(x=gamma,y=loglik))+geom_point()+geom_smooth()

gamma_prof %>% 
  subset(loglik==max(loglik)) %>% unlist() -> mle
simulate(sample,params=mle,nsim=10,as.data.frame=TRUE,include.data=TRUE) -> gammaSim
gammaSim$CumulativeCases[1:39] <- gammaSim$actualCumulative[1:39]
ggplot(data=gammaSim,mapping=aes(x=time,y=CumulativeCases,group=sim,alpha=(sim=="data")))+
scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                    labels=c(`FALSE`="simulation",`TRUE`="data"))+
geom_line()

#Particle Markov Chain Monte Carlo using estimates obtained from mif2
bake(file="CumData-pmcmc.rds",{
  beta_prof %>% ddply(~Beta,subset,loglik==max(loglik)) %>%
    subset(beta_sd > 0.01 & beta_sd < 1.0 & Beta > 1.0 & Beta < 3.0 & mu < 1.0 &
             gamma < 2.0 & rho > 0.05 & rho < 0.99,
           select=-c(loglik,loglik.se)) -> starts
  foreach (start=iter(starts,"row"),.combine=c,
           .options.mpi=list(seed=23781975),
           .packages=c("pomp","magrittr"),.errorhandling="remove") %dopar% 
           {
             sample %>%
               pmcmc(Nmcmc=2000,Np=200,start=unlist(start),
                     proposal=mvn.rw.adaptive(rw.sd=c(Beta=0.3,beta_sd=0.05,mu=0.05,gamma=0.1,rho=0.15, N=500),
                                              scale.start=100,shape.start=100)) -> chain
             chain %>% pmcmc(Nmcmc=10000,proposal=mvn.rw(covmat(chain)))
           }
}) -> chains

chains %>% conv.rec() -> traces
rejectionRate(traces[,c("Beta","beta_sd","mu", "gamma", "rho", "N")])

autocorr.diag(traces[,c("Beta","beta_sd", "mu", "gamma", "rho", "N")])

traces <- window(traces,thin=50,start=1000)
plot(traces[,"Beta"])
plot(traces[,"beta_sd"])
plot(traces[,"mu"])
plot(traces[,"gamma"])
plot(traces[,"rho"])
plot(traces[,"N"])

gelman.diag(traces[,c("Beta","beta_sd","mu", "gamma", "rho", "N")])

summary(traces[,c("Beta","beta_sd","mu", "gamma", "rho", "N")])

theta <- summary(traces)$quantiles[c("Beta","beta_sd","mu","gamma", "rho", "N"),'50%']
simulate(sample,params=theta,nsim=10,as.data.frame=TRUE,include.data=TRUE) -> posteriorSim
posteriorSim$CumulativeCases[1:39] <- posteriorSim$actualCumulative[1:39]
  ggplot(data=posteriorSim,mapping=aes(x=time,y=CumulativeCases,group=sim,alpha=(sim=="data")))+
  scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                     labels=c(`FALSE`="simulation",`TRUE`="data"))+
  geom_line()
  
simulate(sample, params= coef(chains[[1]]), nsim=20,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(mapping=aes(x=time,y=CumulativeCases,group=sim,alpha=(sim=="data")))+
  scale_alpha_manual(name="",values=c(`TRUE`=1,`FALSE`=0.2),
                     labels=c(`FALSE`="simulation",`TRUE`="data"))+
  geom_line()
  
chains %>% filter.traj() %>% melt() %>% 
  subset(rep > 1000 & rep %% 50 == 0) %>%
  dcast(L1+rep+time~variable) %>%
  ddply(~time,summarize,
        prob=c(0.025,0.5,0.975),
        q=quantile(CumulativeCases,prob)) %>% 
  mutate(qq=mapvalues(prob,from=c(0.025,0.5,0.975),to=c("lo","med","hi"))) %>%
  dcast(time~qq,value.var='q') %>% 
  ggplot()+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi),alpha=0.5,fill='blue')+
  geom_line(aes(x=time,y=med),color='blue')+
  geom_point(data=cumulativeSampleData,aes(x=time,y=actualCumulative),color='black',size=2)+
  labs(y="Cases")