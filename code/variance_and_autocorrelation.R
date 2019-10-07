


library(reader)
library(stats)

#load data

#example_cb2 <- read_csv("C:/Users/jjjai/Dropbox (ASU)/NIMBioS_Transient_workgroup/nimbios-transients/data/example_cb2.csv")
#example_jae2 <- read_csv("C:/Users/jjjai/Dropbox (ASU)/NIMBioS_Transient_workgroup/nimbios-transients/data/example_jae2.csv")
#reps <- read_csv("C:/Users/jjjai/Dropbox (ASU)/NIMBioS_Transient_workgroup/nimbios-transients/data/reps.csv")

example_cb2 <- read_csv("data//example_cb2.csv")
example_jae2 <- read_csv("data/example_jae2.csv")
reps <- read_csv("data/reps.csv")

reps_data<-reps[1:100000,2]

reps_data <- t(matrix(reps_data, nrow = 100, ncol = 1000))

##################variance


## for example_cb2
Windows1 <- 10
smooth_example_cb2 <- smth(example_cb2[1:125,2],window = 20,method = "gaussian")
smooth_example_cb2 <- smooth_example_cb2[11:125]
Diff_example_cb2 <- (example_cb2[11:125,2]-smooth_example_cb2)
siz_example_cb2 <- length(Diff_example_cb2)
ending <- siz_example_cb2[1]-Windows1

var_example_cb2 <- matrix(0,ending,1)
for(i in 1:ending){
  var_example_cb2[i,1] <- apply(t(Diff_example_cb2[i:(i+Windows1)]),1,var)  
}

plot(250*var_example_cb2,type="o", pch=22, lty=2, col="blue", ylim=c(0,1.5))
lines(example_cb2[11:(125-Windows1),2],type="o", lty=2, col="red")




## for reps

Windows1 <- 200
Windows2 <- 100
Step <- 10
var_reps <- matrix(0,floor((1000-Windows1-Windows2)/Step),100)
smooth_data <- matrix(0,1000,100)
Diff_data <- matrix(0, 1000-Windows1, 100)
siz_reps <- dim(reps_data)
for(i in 1:siz_reps[2]){
  data_now <- reps_data[1:1000,i]
  smooth_data_now <- smth(data_now,window = Windows1,method = "gaussian")
  smooth_data[1:1000,i] <- smooth_data_now
  smooth_data_now <- smooth_data_now[(Windows1+1):1000]
  Diff_data_now <- (data_now[(Windows1+1):1000]-smooth_data_now)
  Diff_data[1:(1000-Windows1),i] <- Diff_data_now
  siz_data_now <- length(Diff_data_now)
  ending <- floor((siz_data_now[1]-Windows2)/Step)
  for(j in 1:ending){
    var_reps[j,i] <- apply(t(Diff_data_now[((j-1)*Step+1):(j*Step+Windows2)]),1,var)
  }
}

matplot(var_reps, type = "l")

x <- 1:floor((1000-Windows1-Windows2)/Step)
x <- x-1

n <- 68
plot(reps_data[(Windows1+Windows2+1):1000,n],type="o", lty=2, col="red", ylim=c(0.3,1.5))
lines(x*Step, 150*var_reps[1:floor((1000-Windows1-Windows2)/Step),n],type="o", pch=22, lty=2, col="blue")
lines(smooth_data[(Windows1+Windows2+1):1000,n],type="o", lty=2, col="green")
plot(Diff_data_now,type="o", lty=2)


##################Autocorrelation


## for reps

Windows1 <- 200
Windows2 <- 100
Step <- 10
AR_reps <- matrix(0,floor((1000-Windows1-Windows2)/Step),100)
smooth_data <- matrix(0,1000,100)
Diff_data <- matrix(0, 1000-Windows1, 100)
siz_reps <- dim(reps_data)
for(i in 1:siz_reps[2]){
  data_now <- reps_data[1:1000,i]
  smooth_data_now <- smth(data_now,window = Windows1,method = "gaussian")
  smooth_data[1:1000,i] <- smooth_data_now
  smooth_data_now <- smooth_data_now[(Windows1+1):1000]
  Diff_data_now <- (data_now[(Windows1+1):1000]-smooth_data_now)
  Diff_data[1:(1000-Windows1),i] <- Diff_data_now
  siz_data_now <- length(Diff_data_now)
  ending <- floor((siz_data_now[1]-Windows2)/Step)
  for(j in 1:ending){
    model <- Arima(Diff_data_now[((j-1)*Step+1):(j*Step+Windows2)], order = c(0, 0, 0))
    AR_reps[j,i] <- model$coef
  }
}

matplot(AR_reps, type = "l")

x <- 1:floor((1000-Windows1-Windows2)/Step)
x <- x-1

n <- 28
plot(reps_data[(Windows1+Windows2+1):1000,n],type="o", lty=2, col="red", ylim=c(0,1.5))
lines(x*Step, 50*AR_reps[1:floor((1000-Windows1-Windows2)/Step),n],type="o", pch=22, lty=2, col="blue")
lines(smooth_data[(Windows1+Windows2+1):1000,n],type="o", lty=2, col="green")
#plot(Diff_data_now,type="o", lty=2)