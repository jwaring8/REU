flu <- read.csv('flu.csv')
flu
rm(list=ls(all=TRUE))
load('data.RData')
ls()
plot(flu$day,flu$flu,type='b',xlab='Day',ylab='Number of individuals infected')
prevalence<-flu$flu/764
prevalence
plot(flu$day, flu$flu, type='b', xlab='Day', ylab='Number of individuals infected')
names(niamey)[1]<-"District 1"
names(niamey)[2]<-"District 2"
names(niamey)[3]<-"District 3"
View(niamey)
png("Rplot.png")
plot(1:16, niamey$`District 1`, main="Measles Infections in Districts of Niger", type='b', col='green', xlab='Day', ylab='Number of infected')
lines(x=1:16, y = niamey$`District 2`, type = 'b', col='blue')
lines(x=1:16, y = niamey$`District 3`, type='b', col='red')
dev.off() #saving plot as an image called Rplot.png
which(rowSums(niamey, na.rm=T)==max(rowSums(niamey, na.rm =T)))
rowSums(niamey)[10]
biweek.sum<-rowSums(niamey)
biweek<-1:16
plot(biweek, biweek.sum, type='l', ylab='Measles cases in all districts', xlab='Biweek period')
max<-which(biweek.sum==max(biweek.sum, na.rm=T))
abline(v=10, col='red')
