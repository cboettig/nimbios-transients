#' Plot population trajectories
#'
#' @param obs_y observed population with error
#' @param true_x true population
#' @param t.step time step
#'
#' @return NULL
#' @export
plot_traj <- function(obs_y = NULL, true_x = NULL, t.step = 1){
  if(is.null(obs_y)){
    matplot(seq(1, N, l = N/constants$t.step), true_x,
            type = "l", ylab = "x", xlab = "time",
            col = scales::alpha("black", 1), lty = 1)
  } else {
    matplot(seq(1, N, l = N/constants$t.step), obs_y,
            type = "l", ylab = "x", xlab = "time",
            col = scales::alpha("black", 0.4), lty = 1)
    matplot(seq(1, N, l = N/constants$t.step), true_x, type = "l",
            col = scales::alpha("darkred", 0.6), lty = 1, add = T)
    legend("bottomright", col = c("black", "darkred"), lty = 1, legend = c("y", "x"))
  }
  return(NULL)
}
