## Started 30 May 2019 ##
## By Lizzie ##

## Last udpated 11 Septemer 2019 ##

## Next steps are ...
# (1) Check that model is properly coded
# (2) Scale up to test on many runs of data...

library(rstan)
library(nimble)

setwd("~/Documents/git/projects/misc/miscmisc/nimBios2019/nimbios-transients/code/ghost_stan")

if(FALSE){
d <- read.csv("data/example_cb2.csv")
goo <- list(x=d$x, N = nrow(d))
}

# Set up the data
p <- list(r = .05, K = 2, Q = 5, H = .38, sigma = .01, a = 0.023, N = 1e4)
growth <- function(x, p) x * p$r * (1 - x / p$K)
consumption <- function(x,p) p$a * x ^ p$Q / (x^p$Q + p$H^p$Q)
# Define stochastic model in BUGS notation
may  <- nimble::nimbleCode({
  x[1] <- x0
  for(t in 1:(N-1)){
    # Determinstic mean looks like standard R
    mu[t] <- x[t] + x[t] * r * (1 - x[t] / K)  - a * x[t] ^ Q / (x[t] ^ Q + H ^ Q)
    # Note the use of ~ in BUGS to show 'distributed as normal' 
    y[t+1] ~ dnorm(mu[t], sd = sigma)
    x[t+1] <- max(y[t+1],0)
  }
  
})
model <- nimbleModel(may,constants = p, inits = list(x0 = 0.2))
cmodel <- compileNimble(model)
set.seed(123)
simulate(cmodel)

x=cmodel$x
goo <- list(x=x, N = length(x))

goober = stan('ghostfit.stan', data = goo,
               chains=4, iter = 2000, cores=4)

library(shinystan)
launch_shinystan(goober)

sumer <- summary(goober)$summary


# p <- list(r = .05, K = 2, Q = 5, H = .38, sigma = .01, a = 0.023, N = 1e4)

# look at priors
mean(rlnorm(1000, 5, 0.1)) # Q
mean(rlnorm(1000, 0.38, 0.1)) # H
mean(rlnorm(1000, 0.233, 0.1)) # a 
mean(rlnorm(1000, 0.01, 0.001)) # sigma
mean(rlnorm(1000, 0, 1)) # x0

# extra
x2 = x[seq(1, length(x), 5)]
goosmall <- list(x=x2, N = length(x2))

# try to just fit a, fix everything else
goobera = stan('ghostfit_asimple.stan', data = goo,
               chains=4, iter = 1500, cores=4)

sumer_asimple <- summary(goobera)$summary
