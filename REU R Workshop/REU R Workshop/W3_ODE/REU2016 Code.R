#########################
#In class exercise 1: Writing Functions
#########################

ex1 <- function(x){return(5*x^3+6*x-4/x)}

ex1(3)

ex1.x <- seq(-10,-1,0.01)
ex1.y <- ex1(ex1.x)

plot(ex1.x, ex1.y, type="l", xlab="x", ylab="y", lwd=3.0)

ex1.x[which(ex1.y==-7)]

#########################
#In class exercise 2: Discrete vs Continuous Growth
#########################

      par(mfrow=c(3,1), cex.lab=2.0, cex.axis=2.0, cex.main=2.0, mar=c(5,6,4,2))
      geometric.growth1 <- function(N,lambda) {return(N*lambda)}
      
      geometric.growth2 <- function(N0,lambda, t) {return(N0*lambda^t)}
      
      tmax <- 20
      dt <- 1
      t.vec <- seq(0,20,1)
      
      N0 <- 1
      
      N.geo.vec1 <- rep(0,tmax)
      N0 <- 1
      N.geo.vec1[1] <- N0
      
      lambda <- 0.75
      
      for(i in 2:(tmax+1)){
        N.geo.vec1[i] <- geometric.growth1(N.geo.vec1[i-1], lambda)
      }
      
      plot(t.vec, N.geo.vec1, xlab="", ylab="", pch=21, bg="grey40", col="black", cex=2.0, main=expression(lambda ~ "= 0.75"))
      
      lambda <- 1.0
      
      for(i in 2:(tmax+1)){
        N.geo.vec1[i] <- geometric.growth1(N.geo.vec1[i-1], lambda)
      }
      plot(t.vec, N.geo.vec1, xlab="", ylab="N", pch=21, bg="dodgerblue", col="black", cex=2.0, main=expression(lambda ~ "= 1.0"))
      
      lambda <- 2.0
      
      for(i in 2:(tmax+1)){
        N.geo.vec1[i] <- geometric.growth1(N.geo.vec1[i-1], lambda)
      }
      
      plot(t.vec, N.geo.vec1, xlab="time", ylab="", pch=21, bg="firebrick", col="black", cex=2.0, main=expression(lambda ~ "= 2.0"))


      geometric.growth2 <- function(N0,lambda, t) {return(N0*lambda^t)}
      
      par(mfrow=c(3,1), cex.lab=2.0, cex.axis=2.0, cex.main=2.0, mar=c(5,6,4,2))
      tmax <- 20
      dt <- 1
      t.vec <- seq(0,20,1)
      
      N0 <- 1
      
      lambda <- 0.75
      
      N.geo.vec2 <- geometric.growth2(N0, lambda, t.vec)
      
      plot(t.vec, N.geo.vec2, xlab="", ylab="", pch=21, bg="grey40", col="black", cex=2.0, main=expression(lambda ~ "= 0.75"))
      
      lambda <- 1.0
      
      N.geo.vec2 <- geometric.growth2(N0, lambda, t.vec)
      plot(t.vec, N.geo.vec2, xlab="", ylab="N", pch=21, bg="dodgerblue", col="black", cex=2.0, main=expression(lambda ~ "= 1.0"))
      
      lambda <- 2.0
      
      N.geo.vec2 <- geometric.growth2(N0, lambda, t.vec)
      
      plot(t.vec, N.geo.vec2, xlab="time", ylab="", pch=21, bg="firebrick", col="black", cex=2.0, main=expression(lambda ~ "= 2.0"))
      
      
      
      exp.growth.dis <- function(N0, r, t){return(N0*e^(r*t))}
      
      tmax <- 20
      dt <- 1
      t.vec <- seq(0,20,1)
      
      N0 <- 1
      
      r <- -0.05
      
      N.dis.vec <- exp.growth.dis(N0, r, t.vec)
      
      plot(t.vec, N.dis.vec, xlab="", ylab="", pch=21, bg="grey40", col="black", cex=2.0, main=expression(r ~ "= -0.05"))

      r <- 0
      
      N.dis.vec <- exp.growth.dis(N0, r, t.vec)
      
      plot(t.vec, N.dis.vec, xlab="", ylab="N", pch=21, bg="dodgerblue", col="black", cex=2.0, main=expression(r ~ "= 0"))
      
      r <- 0.05
      
      N.dis.vec <- exp.growth.dis(N0, r, t.vec)
      
      plot(t.vec, N.dis.vec, xlab="time", ylab="", pch=21, bg="firebrick", col="black", cex=2.0, main=expression(r ~ "= 0.05"))
      
      

#########################
#In class exercise 3: Exponential Growth
#########################
      
      install.packages("deSolve")
      
      require(deSolve)
      
      ExpGrow <- function(t, y, params){
        
        N <- y["N"]
        
        r <- params["r"]
        
        dN <- r*N
        
        return(list(c(dN)))
        
      }
      
      tmax <- 10
      dt <- 0.1
      t.vec <- seq(0, tmax, dt)  
      
      
      N0 <- 1
      y0 <- c(N = N0)
      
      r <- 0.75
      params <- c(r = r)
      
     odeoutput <- lsoda(y=y0, times=t.vec, func=ExpGrow, parms=params)
     
     odeoutput <- as.data.frame(odeoutput)
     
     par(mfrow=c(1,1), cex.axis=2.0, cex.lab=2.0, pch=21, mar=c(5,5,4,2))
     
     plot(N~time, data=odeoutput, bg="dodgerblue", col="black", cex=2.0)
     
     tail(odeoutput)
     
     tmax <- 10
     dt <- 0.1
     t.vec <- seq(0, tmax, dt)  
     
     N0 <- 100
     y0 <- c(N = N0)
     
     r <- -0.5
     params <- c(r = r)
     
     odeoutput <- lsoda(y=y0, times=t.vec, func=ExpGrow, parms=params)
     
     odeoutput <- as.data.frame(odeoutput)
     
     par(cex.axis=2.0, cex.lab=2.0, pch=21, mar=c(5,5,4,2))
     
     plot(N~time, data=odeoutput, bg="dodgerblue", col="black", cex=2.0)
     
 #########################
 #In class example 1: Closed SIR
 #########################  
     
     
     ClosedSIR <- function(t, y, params){
       
       S <- y["S"]; I <- y["I"]; R <- y["R"]; N <- S+I+R
       
       beta <- params["beta"]; gamma <- params["gamma"]
       
       dS <- -beta*I/N*S
       dI <- beta*I/N*S - gamma*I
       dR <- gamma*I
       
       return(list(c(dS, dI, dR)))
       
     }
     
     tmax <- 120
     dt <- 5
     t.vec <- seq(0, tmax, dt)
     
     beta <- 0.3
     gamma <- 1/7
     
     params <- c(beta=beta, gamma=gamma)
     
     N <- 10000
     I0 <- 1/N
     S0 <- 9999/N
     R0 <- 0/N
     
     y0 <- c(S=S0, I=I0, R=R0)
     
     odeoutput <- lsoda(y=y0, times=t.vec, func=ClosedSIR, parms=params)
     
     odeoutput <- as.data.frame(odeoutput)
     
     plot(S~time, data=odeoutput, type="b", bg="darkgoldenrod1", col="black", ylim=c(0,1), ylab="Proportion")
     lines(I~time, data=odeoutput, type="b", bg="firebrick", col="black")
     lines(R~time, data=odeoutput, type="b", bg="dodgerblue", col="black")
     legend("topright", c("S", "I", "R"), pt.bg=c("darkgoldenrod1", "firebrick", "dodgerblue"), pch=21, cex=2.0)
     
     
     
     plot(I~S, data=odeoutput, type="p", bg="grey40")
     
 #########################
 #In class example 2: Reed Frost
 #########################     
     
     alpha <- 0.96
     
     N <- 100
     i <- 1
     s0 <- N-i
     
     t <- 0
     
     out.rf <- c(t, s0, i)
     
     while(i>0){
       
       s <- rbinom(1, s0, alpha^i)
       i <- s0 - s
       t <- t+1
       
       out.rf <- rbind(out.rf, c(t, s, i))
       
       s0 <- s
       
     }
     
     out.rf <- as.data.frame(out.rf)
     names(out.rf) <- c("time", "S", "I")
     
     plot(S~time, data=out.rf, type="b", bg="darkgoldenrod1", col="black", ylab="Number of Individuals", ylim=c(0,N))
     lines(I~time, data=out.rf, type="b", bg="firebrick", col="black")
     
     legend("topright", c("S", "I"), pt.bg=c("darkgoldenrod1", "firebrick"), pch=21, cex=2.0)
     
     