## ----prelims,include=FALSE,cache=FALSE-----------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
)

set.seed(594709947L)
library(ggplot2)
theme_set(theme_bw())
library(grid)
library(plyr)
library(reshape2)
library(foreach)
library(doParallel)
library(pomp)
stopifnot(packageVersion("pomp")>="1.4.5")

## ----sir-construct-------------------------------------------------------
read.table("http://kingaa.github.io/sbied/stochsim/bsflu_data.txt") -> bsflu

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = nearbyint(N)-1;
  I = 1;
  R = 0;
  H = 0;
")

dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
rmeas <- Csnippet("B = rbinom(H,rho);")

pomp(bsflu,times="day",t0=0,
     rprocess=euler.sim(sir_step,delta.t=1/5),
     initializer=sir_init,rmeasure=rmeas,dmeasure=dmeas,
     zeronames="H",statenames=c("H","S","I","R"),
     paramnames=c("Beta","gamma","rho","N")) -> sir

## ----bbs-mc-like-2,results='markup'--------------------------------------
simulate(sir,params=c(Beta=2,gamma=1,rho=0.8,N=2600),
         nsim=10000,states=TRUE) -> x
matplot(time(sir),t(x["H",1:50,]),type='l',lty=1,
        xlab="time",ylab="H",bty='l',col='blue')
lines(time(sir),obs(sir,"B"),lwd=2,col='black')

## ----bbs-mc-like-3,results='markup',cache=T------------------------------
ell <- dmeasure(sir,y=obs(sir),x=x,times=time(sir),
                params=c(Beta=2,gamma=1,rho=0.8,N=2600),log=TRUE)
dim(ell)

## ----bbs-mc-like-4,results='markup'--------------------------------------
ell <- apply(ell,1,sum); summary(exp(ell)); logmeanexp(ell,se=TRUE)

## ----sir-pfilter-1,results='markup',cache=T------------------------------
pf <- pfilter(sir,Np=5000,params=c(Beta=2,gamma=1,rho=0.8,N=2600))
logLik(pf)

## ----sir-pfilter-2,results='markup',cache=T------------------------------
pf <- replicate(10,pfilter(sir,Np=10000,
                           params=c(Beta=2,gamma=1,rho=0.8,N=2600)))
ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

## ----sir-like-slice,cache=TRUE,results='hide'----------------------------
sliceDesign(
  c(Beta=2,gamma=1,rho=0.8,N=2600),
  Beta=rep(seq(from=0.5,to=4,length=40),each=3),
  gamma=rep(seq(from=0.5,to=2,length=40),each=3)) -> p

library(foreach)
library(doParallel)
registerDoParallel()

set.seed(998468235L,kind="L'Ecuyer")

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
  pfilter(sir,params=unlist(theta),Np=5000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> p

## ----sir-like-slice-plot,cache=FALSE,echo=FALSE,results="hide"-----------
library(magrittr)
library(reshape2)
library(ggplot2)
p %>% 
  melt(measure=c("Beta","gamma")) %>%
  subset(variable==slice) %>%
  ggplot(aes(x=value,y=loglik,color=variable))+
  geom_point()+
  facet_grid(~variable,scales="free_x")+
  guides(color=FALSE)+
  labs(x="parameter value",color="")+
  theme_bw()

## ----sir-grid1-----------------------------------------------------------
expand.grid(Beta=seq(from=1,to=4,length=50),
            gamma=seq(from=0.7,to=3,length=50),
            rho=0.8,
            N=2600) -> p

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
  pfilter(sir,params=unlist(theta),Np=5000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> p

## ----sir-grid1-plot,echo=F,cache=F,purl=T--------------------------------
library(magrittr)
library(reshape2)
library(plyr)
p %>% 
  mutate(loglik=ifelse(loglik>max(loglik)-100,loglik,NA)) %>%
  ggplot(aes(x=Beta,y=gamma,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  geom_contour(color='black',binwidth=3)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(gamma))

## ----sir-partrans--------------------------------------------------------
toEst <- Csnippet("
 TBeta = log(Beta);
 Tgamma = log(gamma);
 Trho = logit(rho);
")

fromEst <- Csnippet("
 TBeta = exp(Beta);
 Tgamma = exp(gamma);
 Trho = expit(rho);
")

pomp(sir,toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     paramnames=c("Beta","gamma","rho")) -> sir

## ----sir-ref-params------------------------------------------------------
coef(sir) <- c(Beta=2,gamma=1,rho=0.8,N=2600)

## ----sir-like-optim-1,echo=T,eval=T,results='markup',cache=T-------------
neg.ll <- function (par, est) {
  allpars <- coef(sir,transform=TRUE)
  allpars[est] <- par
  try(
    freeze(
      pfilter(sir,params=partrans(sir,allpars,dir="fromEst"),
              Np=1000),
      seed=5859684
    )
  ) -> pf
  if (inherits(pf,"try-error")) {
    1e10                 ## a big, bad number
  } else {
    -logLik(pf)
  }
}

## ----sir-like-optim-2,results='markup',cache=T---------------------------
## use Nelder-Mead with fixed RNG seed
fit <- optim(
  par=c(log(1), log(2), log(0.8)),
  est=c("gamma","Beta","rho"),
  fn=neg.ll,
  method="Nelder-Mead",
  control=list(maxit=400,trace=0)
)

mle <- sir
coef(mle,c("gamma","Beta","rho"),transform=TRUE) <- fit$par
coef(mle,c("gamma","Beta","rho")) ## point estimate

fit$val

lls <- replicate(n=5,logLik(pfilter(mle,Np=5000)))
ll <- logmeanexp(lls,se=TRUE); ll

simulate(mle,nsim=8,as.data.frame=TRUE,include.data=TRUE) -> sims

## ----sir-like-optim-plot,echo=F------------------------------------------
ggplot(data=sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  guides(color=FALSE)+
  geom_line()

