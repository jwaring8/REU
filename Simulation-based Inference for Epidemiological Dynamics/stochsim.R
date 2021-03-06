## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
  )

## ----prelims,echo=F,cache=F----------------------------------------------
set.seed(594709947L)
library(plyr)
library(reshape2)
library(pomp)
library(ggplot2)
theme_set(theme_bw())
stopifnot(packageVersion("pomp")>="1.4.5")

## ----flu-data1-----------------------------------------------------------
read.table("http://kingaa.github.io/sbied/stochsim/bsflu_data.txt") -> bsflu
head(bsflu)

## ----flu-data2,echo=F----------------------------------------------------
ggplot(data=bsflu,aes(x=day,y=B))+geom_line()+geom_point()

## ----rproc1--------------------------------------------------------------
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
")

sir_init <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
")

pomp(bsflu,time="day",t0=0,rprocess=euler.sim(sir_step,delta.t=1/6),
     initializer=sir_init,paramnames=c("N","Beta","gamma"),
     statenames=c("S","I","R")) -> sir

## ----rproc2--------------------------------------------------------------
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  H = 0;
")

pomp(sir,rprocess=euler.sim(sir_step,delta.t=1/6),initializer=sir_init,
     paramnames=c("Beta","gamma","N"),statenames=c("S","I","R","H")) -> sir

## ----zero1---------------------------------------------------------------
pomp(sir,zeronames="H") -> sir

## ----meas-model----------------------------------------------------------
dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
rmeas <- Csnippet("B = rbinom(H,rho);")

## ----add-meas-model------------------------------------------------------
sir <- pomp(sir,rmeasure=rmeas,dmeasure=dmeas,statenames="H",paramnames="rho")

## ------------------------------------------------------------------------
sims <- simulate(sir,params=c(Beta=1.5,gamma=1,rho=0.9,N=2600),
                 nsim=20,as=TRUE,include=TRUE)

ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

