require(ggplot2)
require(gridExtra)
require(plyr)

load("~/R/SimulatedData/SimulatedData/allRawTrajMatches.Rda")
load("~/R/SimulatedData/SimulatedData/simulations.Rda")
load("~/R/SimulatedData/SimulatedData/allCumTrajMatches.Rda")

#Getting data to plot all 25 raw simulations together
simPlots <- list()
for(i in 1:24){
  sampleData <- allSims[[i]]
  plot <- ggplot(data=sampleData, mapping=aes(x=time, y= Cases)) +
    geom_line()+expand_limits(y=0)+theme_classic()
  simPlots[[length(simPlots)+1]] <- plot 
}
setwd("~/R/SimulatedData/SimulatedData")
sampleData <- read.csv("sample_data.csv")
myplot1<-ggplot(data=sampleData, mapping=aes(x=time, y= Cases)) +
  geom_line()+expand_limits(y=0)+theme_classic()
myplot2<-simPlots[[1]]
myplot3<-simPlots[[2]]
myplot4<-simPlots[[3]]
myplot5<-simPlots[[4]]
myplot6<-simPlots[[5]]
myplot7<-simPlots[[6]]
myplot8<-simPlots[[7]]
myplot9<-simPlots[[8]]
myplot10<-simPlots[[9]]
myplot11<-simPlots[[10]]
myplot12<-simPlots[[11]]
myplot13<-simPlots[[12]]
myplot14<-simPlots[[13]]
myplot15<-simPlots[[14]]
myplot16<-simPlots[[15]]
myplot17<-simPlots[[16]]
myplot18<-simPlots[[17]]
myplot19<-simPlots[[18]]
myplot20<-simPlots[[19]]
myplot21<-simPlots[[20]]
myplot22<-simPlots[[21]]
myplot23<-simPlots[[22]]
myplot24<-simPlots[[23]]
myplot25<-simPlots[[24]]

grid.arrange(myplot1, myplot2, myplot3, myplot4, myplot5, myplot6, myplot7, myplot8, myplot9, myplot10,
             myplot11, myplot12, myplot13, myplot14, myplot15, myplot16, myplot17, myplot18, myplot19, myplot20,
             myplot21, myplot22, myplot23, myplot24, myplot25, ncol=5)

#Getting data to plot all 25 cumulative simulations together
sampleData <- read.csv(file ="sample_data.csv")
sampleData <- subset(sampleData, select=c("time", "Cases"))
gamma <- 0.6
cumulativeSampleData <- mutate(sampleData, Cases=trunc(cumsum(Cases-(gamma*Cases))))
myplot1<-ggplot(data=cumulativeSampleData, mapping=aes(x=time, y= Cases)) +
  geom_line()+expand_limits(y=0)+theme_classic()
simPlots <- list()
for(i in 1:24){
  sampleData <- allSims[[i]]
  sampleData <- subset(sampleData, select=c("time", "Cases"))
  cumulativeSampleData <- mutate(sampleData, Cases=trunc(cumsum(Cases-(gamma*Cases))))
  ggplot(data=cumulativeSampleData, mapping=aes(x=time, y= Cases)) +
    geom_line()+expand_limits(y=0)+theme_classic() -> plot
  simPlots[[length(simPlots)+1]] <- plot 
}
myplot2<-simPlots[[1]]
myplot3<-simPlots[[2]]
myplot4<-simPlots[[3]]
myplot5<-simPlots[[4]]
myplot6<-simPlots[[5]]
myplot7<-simPlots[[6]]
myplot8<-simPlots[[7]]
myplot9<-simPlots[[8]]
myplot10<-simPlots[[9]]
myplot11<-simPlots[[10]]
myplot12<-simPlots[[11]]
myplot13<-simPlots[[12]]
myplot14<-simPlots[[13]]
myplot15<-simPlots[[14]]
myplot16<-simPlots[[15]]
myplot17<-simPlots[[16]]
myplot18<-simPlots[[17]]
myplot19<-simPlots[[18]]
myplot20<-simPlots[[19]]
myplot21<-simPlots[[20]]
myplot22<-simPlots[[21]]
myplot23<-simPlots[[22]]
myplot24<-simPlots[[23]]
myplot25<-simPlots[[24]]

grid.arrange(myplot1, myplot2, myplot3, myplot4, myplot5, myplot6, myplot7, myplot8, myplot9, myplot10,
             myplot11, myplot12, myplot13, myplot14, myplot15, myplot16, myplot17, myplot18, myplot19, myplot20,
             myplot21, myplot22, myplot23, myplot24, myplot25, ncol=5)

#Plot all 25 simulations for each of the 38 traj matches together (RAW)
for(i in 1:38) {
  myplot1<-allRawTrajMatches[[1]][[i]]
  myplot2<-allRawTrajMatches[[2]][[i]]
  myplot3<-allRawTrajMatches[[3]][[i]]
  myplot4<-allRawTrajMatches[[4]][[i]]
  myplot5<-allRawTrajMatches[[5]][[i]]
  myplot6<-allRawTrajMatches[[6]][[i]]
  myplot7<-allRawTrajMatches[[7]][[i]]
  myplot8<-allRawTrajMatches[[8]][[i]]
  myplot9<-allRawTrajMatches[[9]][[i]]
  myplot10<-allRawTrajMatches[[10]][[i]]
  myplot11<-allRawTrajMatches[[11]][[i]]
  myplot12<-allRawTrajMatches[[12]][[i]]
  myplot13<-allRawTrajMatches[[13]][[i]]
  myplot14<-allRawTrajMatches[[14]][[i]]
  myplot15<-allRawTrajMatches[[15]][[i]]
  myplot16<-allRawTrajMatches[[16]][[i]]
  myplot17<-allRawTrajMatches[[17]][[i]]
  myplot18<-allRawTrajMatches[[18]][[i]]
  myplot19<-allRawTrajMatches[[19]][[i]]
  myplot20<-allRawTrajMatches[[20]][[i]]
  myplot21<-allRawTrajMatches[[21]][[i]]
  myplot22<-allRawTrajMatches[[22]][[i]]
  myplot23<-allRawTrajMatches[[23]][[i]]
  myplot24<-allRawTrajMatches[[24]][[i]]
  myplot25<-allRawTrajMatches[[25]][[i]]
  
  grid.arrange(myplot1, myplot2, myplot3, myplot4, myplot5, myplot6, myplot7, myplot8, myplot9, myplot10,
               myplot11, myplot12, myplot13, myplot14, myplot15, myplot16, myplot17, myplot18, myplot19, myplot20,
               myplot21, myplot22, myplot23, myplot24, myplot25, ncol=5)
}

#Plot all 25 simulations for each of the 38 traj matches together (CUM)
for(i in 2:38) {
  myplot1<-allCumTrajMatches[[1]][[i]]
  myplot2<-allCumTrajMatches[[2]][[i]]
  myplot3<-allCumTrajMatches[[3]][[i]]
  myplot4<-allCumTrajMatches[[4]][[i]]
  myplot5<-allCumTrajMatches[[5]][[i]]
  myplot6<-allCumTrajMatches[[6]][[i]]
  myplot7<-allCumTrajMatches[[7]][[i]]
  myplot8<-allCumTrajMatches[[8]][[i]]
  myplot9<-allCumTrajMatches[[9]][[i]]
  myplot10<-allCumTrajMatches[[10]][[i]]
  myplot11<-allCumTrajMatches[[11]][[i]]
  myplot12<-allCumTrajMatches[[12]][[i]]
  myplot13<-allCumTrajMatches[[13]][[i]]
  myplot14<-allCumTrajMatches[[14]][[i]]
  myplot15<-allCumTrajMatches[[15]][[i]]
  myplot16<-allCumTrajMatches[[16]][[i]]
  myplot17<-allCumTrajMatches[[17]][[i]]
  myplot18<-allCumTrajMatches[[18]][[i]]
  myplot19<-allCumTrajMatches[[19]][[i]]
  myplot20<-allCumTrajMatches[[20]][[i]]
  myplot21<-allCumTrajMatches[[21]][[i]]
  myplot22<-allCumTrajMatches[[22]][[i]]
  myplot23<-allCumTrajMatches[[23]][[i]]
  myplot24<-allCumTrajMatches[[24]][[i]]
  myplot25<-allCumTrajMatches[[25]][[i]]
  
  grid.arrange(myplot1, myplot2, myplot3, myplot4, myplot5, myplot6, myplot7, myplot8, myplot9, myplot10,
               myplot11, myplot12, myplot13, myplot14, myplot15, myplot16, myplot17, myplot18, myplot19, myplot20,
               myplot21, myplot22, myplot23, myplot24, myplot25, ncol=5)
}

