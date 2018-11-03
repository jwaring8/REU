#simple stochastic SIR model with births and deaths


SIR_simulation <-function (gamma,mu,beta,N,dt,number.timeunits){
	model.output=matrix(0,nrow=3,ncol=number.timeunits)

	#initialise model
	#row1 = S
	#row2 = I
	#row3 = R

	model.output[1,1]=N-1
	model.output[2,1]=1 #assume one infected individual at start

	#run stochastic model
	for(i in 2:number.timeunits){
		#evaluate random number leaving S(incidence and death)
		prob.inc = dt*beta*model.output[2,(i-1)]/(sum(model.output[,(i-1)]))
		prob.death = mu*dt
		inc.death.remainS = rmultinom(1,model.output[1,(i-1)],c(
				prob.inc,prob.death,1-(prob.inc + prob.death)))

		#evaluate random number leaving I (recovery and death)
		prob.rec = dt * gamma
	
		rec.death.remainI = rmultinom(1,model.output[2,(i-1)],c(prob.rec,prob.death,1-(prob.rec + prob.death)))

		#evaluate random number leaving R
		death.remainR = rmultinom(1,model.output[3,(i-1)],c(prob.death,1-prob.death))

		#evaluate number of births
		births = rbinom(1,sum(model.output[,(i-1)]),mu*dt)

		#update timestep
		model.output[1,i] = model.output[1,(i-1)] + births - sum(inc.death.remainS[1:2])
		model.output[2,i] = model.output[2,(i-1)] + inc.death.remainS[1] - sum(rec.death.remainI[1:2])
		model.output[3,i] = model.output[3,(i-1)] + rec.death.remainI[1] - death.remainR[1]
  

	}
	return(model.output)
}

R0=1.2
gamma = 1/3 #(3 days duration infection)

mu = 1/(60*365) #(life expectancy = 60 years)

beta = R0/(gamma + mu)

N=1000

dt = 0.01

number.timeunits = 365*0.1*(1/dt)

onesim=SIR_simulation(gamma,mu,beta,N,dt,number.timeunits)

plot((1:number.timeunits)*dt,colSums(onesim),type="l",ylim=c(0,N*1.1)) #total population - should be approximately constant
lines((1:number.timeunits)*dt,onesim[1,],col="blue")#S
lines((1:number.timeunits)*dt,onesim[2,],col="red")#I
lines((1:number.timeunits)*dt,onesim[3,],col="purple")#R

#100 simulations
nsims=100
results=array(0,c(3,number.timeunits,nsims))

for(i in 1:nsims){
	results[,,i]=SIR_simulation(gamma,mu,beta,N,dt,number.timeunits)
}
plot((1:number.timeunits)*dt,results[2,,1],type="l",ylim=c(0,N*1.1),xlab="days",ylab="Number infected") #total population - should be approximately constant
for(i in 2:nsims){
	lines((1:number.timeunits)*dt,results[2,,i]) #total population - should be approximately constant
}