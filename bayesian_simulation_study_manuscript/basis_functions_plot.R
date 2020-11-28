## basis functions
x <- seq(0.2, 1.8, l = 2e2)
X <- cbind(1.524111,
           7.180e-02 * x - 7.495e-02,
           1.823e-01 * x^2 - 3.805e-01 * x + 1.633e-01,
           4.712e-01 * x^3 - 1.476e+00 * x^2 + 1.376e+00 * x - 3.642e-01,
           1.226 * x^4 - 5.118 * x^3 + 7.402 * x^2 - 4.300 * x + 8.247e-01)
pdf("fig/basis_functions.pdf", height = 6)
matplot(x, X[, -1], type = "l", col = "gray", lty = 2,
        xlab = "population", ylab = expression(phi(x)), main = "polynomial basis functions")
dev.off()