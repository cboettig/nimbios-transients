## libraries ----
library(parallel)
library(nimble)
## parameters ----
N <- 1e3
r <- 0.05; K <- 2
a <- 0.023; H <- 0.38; Q <- 5
x0 <- 0.3; sigma <- 0.02
## define truncated normal distribution ----
dtruncnorm <- nimbleFunction(
  run = function(x = double(0), mean = double(0), 
                 sd = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    log_prob <- dnorm(x = x, mean = mean, sd = sd, log = 1) - 
      pnorm(q = 0, mean = -mean, sd = sd, log = 1)
    if(log) return(log_prob)
    else return(exp(log_prob))
  })
rtruncnorm <- nimbleFunction(
  run = function(n = integer(0, default = 1), mean = double(0),
                 sd = double(0)) {
    returnType(double(0))
    if(n != 1) print("rtruncnorm only allows n = 1; using n = 1.")
    draw <- rnorm(n = 1, mean = mean, sd = sd)
    while(draw < 0)
      draw <- rnorm(n = 1, mean = mean, sd = sd)
    return(draw)
  })
## define stochastic model in BUGS notation ----
code <- nimble::nimbleCode({
  log(r) ~ dnorm(mu_r, sd_r)
  log(K) ~ dnorm(mu_K, sd_K)
  log(a) ~ dnorm(mu_a, sd_a)
  log(H) ~ dnorm(mu_H, sd_H)
  log(Q) ~ dnorm(mu_Q, sd_Q)
  log(sigma) ~ dnorm(mu_sigma, sd_sigma)
  x[1] <- x0
  for(t in 1:((1/t.step)*N-1)){
    # Determinstic mean looks like standard R
    mu[t] <- x[t] + t.step*(x[t] * r * (1 - x[t] / K)  - a * x[t] ^ Q / (x[t] ^ Q + H ^ Q))
    # Note the use of ~ in BUGS to show 'distributed as normal' 
    sd_y[t] <- sigma*mu[t]*sqrt(t.step)
    y[t+1] ~ dnorm(mu[t], sd = sd_y[t]) # or should this be lognormal?
    # note: since variance scales linearly with time, sd scales with the square root of time
    x[t+1] <- max(y[t+1],0)
  }
})