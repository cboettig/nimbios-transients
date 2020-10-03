library(FNN)
library(xtable)
source("code/scharf_nimble/model.R")

seed <- 270

# posteriors --------------------------------------------------------------
load(paste0("data/scharf_nimble/samples_", 1, "_", seed, ".RData"))
samples_1 <- samples
load(paste0("data/scharf_nimble/samples_", 10, "_", seed, ".RData"))
samples_10 <- samples
load(paste0("data/scharf_nimble/samples_", 100, "_", seed, ".RData"))
samples_100 <- samples
rm(samples)
nmcmc <- nrow(samples_1)

# priors ------------------------------------------------------------------
priors <- samples_1 * NA
priors[, "log_r"] <- rnorm(nmcmc, log(r), 1)
priors[, "log_K"] <- rnorm(nmcmc, log(K), 1)
priors[, "log_a"] <- rnorm(nmcmc, log(a), 1)
priors[, "log_H"] <- rnorm(nmcmc, log(H), 1)
priors[, "log_Q"] <- rnorm(nmcmc, log(Q), 1)
priors[, "log_sigma"] <- rnorm(nmcmc, log(sigma), 1)

# KL ----------------------------------------------------------------------
k <- 10 # KL divergence max number of nearest neighbors to search
par_nm <- c("log_r", "log_K", "log_a", "log_H", "log_Q", "log_sigma")
kl_tab <- matrix(
  nrow = length(par_nm),
  ncol = 3,
  dimnames = list(
    par_nm,
    c("post1", "post10", "post100")
  )
)

for (i in seq_along(par_nm)) {
  kl_tab[par_nm[i], "post1"] <- KL.divergence(
    priors[, par_nm[i]],
    samples_1[, par_nm[i]],
    k
  )[k]
  kl_tab[par_nm[i], "post10"] <- KL.divergence(
    priors[, par_nm[i]],
    samples_10[, par_nm[i]],
    k
  )[k]
  kl_tab[par_nm[i], "post100"] <- KL.divergence(
    priors[, par_nm[i]],
    samples_100[, par_nm[i]],
    k
  )[k]
}

kl_tab

xtable(kl_tab) # create latex table
