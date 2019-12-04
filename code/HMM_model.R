library(depmixS4)
library(ggplot2)

#load observations
#observations <- read.csv(paste(getwd(),"/data/example_jae2.csv",sep=""))

#construct the HMM based on Chen et al 2016
mod <- depmix(x~1,data=observations,nstates=2,trstart=runif(4),instart=c(1,0))
fm <- fit(mod)
results <- fm@posterior
results$state <- as.factor(results$state)
ggplot(data=cbind(observations,results),aes(x=t,y=x,colour=state))+geom_line()

