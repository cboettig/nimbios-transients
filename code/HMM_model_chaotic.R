library(depmixS4)
library(dplyr)
library(ggplot2)

#load observations
observations <- read.csv(paste(getwd(),"/data/HPM_chaotic_transient_data_1.csv",sep=""))
colnames(observations) <- c('t','X1','X2','X3')
observations <- observations %>% filter(t %in% seq(from=1,to=40000,by=1))

#construct the HMM based on Chen et al 2016
mod <- depmix(X1~1,data=observations,nstates=2,trstart=runif(4),instart=c(1,0))
fm <- fit(mod)
results <- fm@posterior

visualizationData <- cbind(observations,results)
ggplot(visualizationData,aes(t,X1))+geom_line(color=as.factor(visualizationData$state))
