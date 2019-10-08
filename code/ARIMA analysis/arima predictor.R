library(forecast)

rawData <-
  read.csv("~/GitHub/nimbios-transients/data/example_jae2.csv")
tsData <- ts(rawData$x[1:15])
fittedData <- auto.arima(tsData)
futureVal <- forecast(fittedData,h=15,levels=c(80,95))

plot(futureVal,ylim=c(0,15))
lines(rawData$x[1:30])