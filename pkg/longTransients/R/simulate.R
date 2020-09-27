#' Simulate population trajectories 
#'
#' @param constants list with constants for nimbleModel()
#' @param inits list with initial values for nimbleModel()
#' @param INCLUDE_ME logical, default \code{TRUE}; should error in population measurements be included
#' @param seed 
#' @importFrom nimble nimbleFunction
#'
#' @return list with population trajectory and observations with error if applicable; also constants and initial values used to simulate 
#' @export
simulate_model_nM <- function(constants = NULL, inits = NULL, INCLUDE_ME = T, seed = 270){
  require(nimble)
  code <- make_model_code_nC()
  model <- nimbleModel(code = code, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  ## set seed + simulate ----
  set.seed(seed)
  simulate(cmodel, nodes = c("y", "x", "mu", "sd_x"))
  true_x <- cmodel$x
  if(INCLUDE_ME) obs_y <- cmodel$y
  return(list(obs_y = obs_y, true_x = true_x, constants = constants, inits = inits))
}
