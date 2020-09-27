.onLoad <- function(libname, longTransients) {
dspikenorm <<- nimbleFunction(
  run = function(x = double(0), mean = double(0), 
                 sd = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if(x < 0){
      if(log) return(-Inf)
      else return(0)
    }
    if(x == 0){
      log_prob <- pnorm(q = 0, mean = mean, sd = sd, log = 1)
      if(log) return(log_prob)
      else return(exp(log_prob))
    }
    log_prob <- dnorm(x = x, mean = mean, sd = sd, log = 1)
    if(log) return(log_prob)
    else return(exp(log_prob))
  })
rspikenorm <<- nimbleFunction(
  run = function(n = integer(0, default = 1), mean = double(0),
                 sd = double(0)) {
    returnType(double(0))
    if(n != 1) print("rtruncnorm only allows n = 1; using n = 1.")
    draw <- rnorm(n = 1, mean = mean, sd = sd)
    if(draw <= 0) return(0)
    else return(draw)
  })
dV <<- nimbleFunction(
  run = function(x = double(0), beta = double(1)) {
    returnType(double(0))
    X <- c(1.524111, 
           7.180e-02 * x - 7.495e-02,
           1.823e-01 * x^2 - 3.805e-01 * x + 1.633e-01,
           4.712e-01 * x^3 - 1.476e+00 * x^2 + 1.376e+00 * x - 3.642e-01,
           1.226 * x^4 - 5.118 * x^3 + 7.402 * x^2 - 4.300 * x + 8.247e-01)
    dV <- inprod(X, beta)
    return(dV)
  })
}