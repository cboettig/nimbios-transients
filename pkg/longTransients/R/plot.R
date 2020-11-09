#' Plot population trajectories
#'
#' @param obs_y observed population with error
#' @param true_x true population
#' @param t.step time step
#'
#' @importFrom scales alpha
#'
#' @return NULL
#' @export
plot_traj <- function(obs_y = NULL, true_x = NULL, t.step = 1){
  N <- nrow(true_x)
  if(is.null(obs_y)){
    matplot(seq(1, N, l = N / t.step), true_x,
            type = "l", ylab = "x", xlab = "time",
            col = alpha("black", 1), lty = 1)
  } else {
    matplot(seq(1, N, l = N / t.step), obs_y,
            type = "l", ylab = "x", xlab = "time",
            col = alpha("black", 0.4), lty = 1)
    matplot(seq(1, N, l = N / t.step), true_x, type = "l",
            col = alpha("darkred", 0.6), lty = 1, add = T)
    legend("bottomright", col = c("black", "darkred"), lty = 1, legend = c("y", "x"))
  }
  return(NULL)
}

#' Traceplots for MCMC output from parameteric model
#'
#' @param samples 
#' @param true values used to generate data for comparison
#'
#' @return NULL
#' @export
plot_trace <- function(samples, true = NULL){
  variables <- c("r", "K", "a", "H", "Q", "sigma")
  INCLUDE_ME <- "sigma_me" %in% colnames(samples)
  if(INCLUDE_ME) variables <- c(variables, "sigma_me")
  layout(matrix(1:8, 4, 2))
  par(mar = c(2, 2, 4, 2))
  for(var in variables){
    plot(samples[, var], type = "l", main = var, ylab = "")
    abline(h = true[var])
  }
  corrplot::corrplot(cor(samples[, variables]))
  return(NULL)
}

#' Traceplots for MCMC output from functional model
#'
#' @param samples 
#'
#' @return NULL
#' @export
plot_trace_functional <- function(samples){
  layout(matrix(c(1:3, 3), 2, 2), widths = c(1, 2))
  variables <- c(paste0("beta[", 1:5, "]"), "sigma", "sigma_me")
  matplot(samples[, paste0("beta[", 1:5, "]")], type = "l")
  matplot(samples[, c("sigma", "sigma_me")], type = "l", lty = 1)
  legend("bottomright", col = 1:2, lty = 1, bty = "n", 
         legend = c(expression(sigma), expression(sigma[me])))
  corrplot::corrplot(cor(samples[, variables]))
  return(NULL)
}

#' Plot realizations of the potential function and its derivative for the parametric model
#'
#' @param samples 
#' @param n_subset 
#' @param true 
#' @param obs_y observations
#' @param x 
#' @param probs 
#' @param alpha 
#' @param ylim_potential 
#' @param ylim_grad 
#' @param which_plot \code{"potential"} for potential function, "grad" for negative gradient of potential function, or \code{"both"} for both
#' @param main 
#'
#' @return
#' @export
plot_potential <- function(samples, n_subset = 400, true = NULL, obs_y = NA,
                           x = seq(0, 2, length.out = 1e2),
                           probs = c(0.125, 0.5, 0.875), alpha = 2e-2,
                           ylim_potential = NULL, ylim_grad = NULL,
                           which_plot = "both", main = NULL){
  subset <- sample(1:nrow(samples), min(n_subset, nrow(samples)))
  dpotential_curves <- apply(samples[subset, ], 1, function(row){
    sapply(x, dpotential, a = row['a'], r = row['r'], 
           H = row['H'], Q = row['Q'], K = row['K'])
  })
  potential_curves <- apply(-dpotential_curves[-1, ] * diff(x), 2, cumsum)
  dpotential_curves_quantiles <- apply(dpotential_curves, 1, quantile, probs = probs)
  potential_curves_quantiles <- apply(potential_curves, 1, quantile, probs = probs)
  true_dpotential <- dpotential(x = x, a = true$a, r = true$r, 
                                H = true$H, Q = true$Q, K = true$K)
  true_potential <- potential(x = x, a = true$a, r = true$r, 
                              H = true$H, Q = true$Q, K = true$K)
  if(which_plot == "both") layout(matrix(1:2, 1, 2))
  if(which_plot %in% c("potential", "both")){
    if(is.null(main)) main <- "potential function"
    matplot(x[-1], t(potential_curves_quantiles), type = "l", 
            lty = c(2, 1, 2), lwd = c(1, 2, 1), 
            main = main, ylim = ylim_potential,
            ylab = "", xlab = "population", col = "black")
    matplot(x[-1], potential_curves, type = "l", lty = 1, 
            col = scales::alpha("black", alpha), lwd = 2, add = T)
    lines(x[-1], true_potential, lwd = 2, col = "darkred")
    rug(obs_y)
  }
  if(which_plot %in% c("grad", "both")){
    if(is.null(main)) main <- "gradient of potential function"
    matplot(x, -t(dpotential_curves_quantiles), type = "l", 
            lty = c(2, 1, 2), lwd = c(1, 2, 1), 
            ylim = ylim_grad, main = main, 
            ylab = "", xlab = "population", col = "black")
    matplot(x, -dpotential_curves, type = "l", lty = 1,
            col = scales::alpha("black", alpha), lwd = 2, add = T)
    lines(x, -true_dpotential, lwd = 2, col = "darkred")
    abline(h = 0, lwd = 2, lty = 3)
    rug(obs_y)
  }
}

#' Plot realizations of the potential function and its derivative for the parametric model
#'
#' @param samples 
#' @param n_subset 
#' @param true 
#' @param obs_y observations
#' @param x 
#' @param probs 
#' @param alpha 
#' @param ylim_potential 
#' @param ylim_grad 
#' @param which_plot \code{"potential"} for potential function, "grad" for negative gradient of potential function, or \code{"both"} for both
#' @param main 
#'
#' @return
#' @export
plot_potential_functional <- function(samples, n_subset = 400, true = NULL, obs_y = NULL,
                                      x = seq(0, 2, length.out = 1e2),
                                      probs = c(0.125, 0.5, 0.875), alpha = 2e-2,
                                      ylim_potential = NULL, ylim_grad = NULL, 
                                      which_plot = "both", main = NULL){
  subset <- sample(1:nrow(samples), min(n_subset, nrow(samples)))
  degree <- length(grep("beta", colnames(samples)))
  dpotential_curves <- apply(samples[subset, 1:degree], 1, function(row){
    sapply(x, dV, beta = row)
  })
  potential_curves <- apply(-dpotential_curves[-1, ] * diff(x), 2, cumsum)
  dpotential_curves_quantiles <- apply(dpotential_curves, 1, quantile, probs = probs)
  potential_curves_quantiles <- apply(potential_curves, 1, quantile, probs = probs)
  true_dpotential <- dpotential(x = x, a = true$a, r = true$r, 
                                H = true$H, Q = true$Q, K = true$K)
  true_potential <- potential(x = x, a = true$a, r = true$r, 
                              H = true$H, Q = true$Q, K = true$K)
  if(which_plot == "both") layout(matrix(1:2, 1, 2))
  if(which_plot %in% c("potential", "both")){
    if(is.null(main)) main <- "potential function"
    matplot(x[-1], t(potential_curves_quantiles), type = "l", 
            lty = c(2, 1, 2), lwd = c(1, 2, 1), 
            main = main, ylim = ylim_potential,
            ylab = "", xlab = "population", col = "black")
    matplot(x[-1], potential_curves, type = "l", lty = 1, 
            col = scales::alpha("black", alpha), lwd = 2, add = T)
    lines(x[-1], true_potential, lwd = 2, col = "darkred")
    rug(obs_y)
  }
  if(which_plot %in% c("grad", "both")){
    if(is.null(main)) main <- "gradient of potential function"
    matplot(x, -t(dpotential_curves_quantiles), type = "l", 
            lty = c(2, 1, 2), lwd = c(1, 2, 1), 
            ylim = ylim_grad, main = main, 
            ylab = "", xlab = "population", col = "black")
    matplot(x, -dpotential_curves, type = "l", lty = 1,
            col = scales::alpha("black", alpha), lwd = 2, add = T)
    lines(x, -true_dpotential, lwd = 2, col = "darkred")
    abline(h = 0, lwd = 2, lty = 3)
    rug(obs_y)
  }
}