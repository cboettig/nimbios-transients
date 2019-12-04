## growth + consumption ----
growth <- function(x, r, K){x * r * (1 - x / K)}
consumption <- function(x, a, H, Q){a * x^Q / (x^Q + H^Q)}
## trace plots ----
layout(matrix(1:8, 4, 2))
par(mar = c(2, 2, 5, 2))
plot(exp(samples[, 'log_r']), type = "l", main = "r", ylab = "")
abline(h = r)
plot(exp(samples[, 'log_K']), type = "l", main = "K", ylab = "")
abline(h = K)
plot(exp(samples[, 'log_a']), type = "l", main = "a", ylab = "",
     ylim = range(exp(samples[, 'log_a']), a))
abline(h = a)
plot(exp(samples[, 'log_H']), type = "l", main = "H", ylab = "")
abline(h = H)
plot(exp(samples[, 'log_Q']), type = "l", main = "Q", ylab = "", 
     log = "y")
abline(h = Q)
plot(exp(samples[, 'log_sigma']), type = "l", main = "sigma", ylab = "")
abline(h = sigma)
plot(samples[, 'mu0'], type = "l", main = "mu0", ylab = "")
abline(h = mu0)
corrplot::corrplot(cor(samples[, c("log_r", "log_K", "log_a", "log_H", "log_Q", "log_sigma", "mu0")]))