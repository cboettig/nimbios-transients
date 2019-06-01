greta.Rmd
================
Carl Boettiger
5/30/2019

``` r
library(tidyverse)
library(greta) # remotes::install_github("greta-dev/greta")
```

``` r
data <- read_csv("../data/reps.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   t = col_double(),
    ##   x = col_double(),
    ##   reps = col_double()
    ## )

``` r
data %>% ggplot(aes(t,x, group=reps)) + geom_line(alpha=0.1)
```

![](greta_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
wide <- data %>% spread(reps, x) %>% select(-t) %>% select(1:5) %>% as.matrix()
n <- dim(wide)[1]
```

``` r
x_t1 <- wide[-1,] # data$x[-1]
x_t <- wide[-n,] # data$x[-100]
dim(x_t)
```

    ## [1] 999   5

``` r
dim(x_t1)
```

    ## [1] 999   5

``` r
#tibble(x_t, x_t1) %>% ggplot(aes(x_t, x_t1)) + geom_point()
```

``` r
#r <- 0.05 
Q <- 5
#sigma <- 0.02
#a <- 0.0233
K <- uniform(0, 10)
H <- uniform(0,  2)
a <- uniform(0, 1)
r <- uniform(0, 1)
#Q <- uniform(0, 10)
sigma <- uniform(0, 1)

# Model
mean <- x_t + r * x_t * (1 - x_t / K) - a * x_t ^ Q / (x_t ^ Q + H ^ Q)
distribution(x_t1) <- normal(mean, sigma)
m <- model(K, a, H, sigma, r)
```

``` r
draws <- mcmc(m, n_samples = 1000, warmup = 3000, chains = 4)
```

    ## 
    ## running 4 chains simultaneously on up to 4 cores

    ## 
        warmup                                           0/3000 | eta:  ?s          
        warmup =                                        50/3000 | eta:  4m | 4% bad 
        warmup =                                       100/3000 | eta:  4m | 2% bad 
        warmup ==                                      150/3000 | eta:  4m | 1% bad 
        warmup ===                                     200/3000 | eta:  3m | 1% bad 
        warmup ===                                     250/3000 | eta:  3m | <1% bad
        warmup ====                                    300/3000 | eta:  3m | <1% bad
        warmup ====                                    350/3000 | eta:  3m | <1% bad
        warmup =====                                   400/3000 | eta:  3m | <1% bad
        warmup ======                                  450/3000 | eta:  3m | <1% bad
        warmup ======                                  500/3000 | eta:  3m | <1% bad
        warmup =======                                 550/3000 | eta:  3m | <1% bad
        warmup ========                                600/3000 | eta:  3m | <1% bad
        warmup ========                                650/3000 | eta:  2m | <1% bad
        warmup =========                               700/3000 | eta:  2m | <1% bad
        warmup ==========                              750/3000 | eta:  2m | <1% bad
        warmup ==========                              800/3000 | eta:  2m | <1% bad
        warmup ===========                             850/3000 | eta:  2m | <1% bad
        warmup ===========                             900/3000 | eta:  2m | <1% bad
        warmup ============                            950/3000 | eta:  2m | <1% bad
        warmup =============                          1000/3000 | eta:  2m | <1% bad
        warmup =============                          1050/3000 | eta:  2m | <1% bad
        warmup ==============                         1100/3000 | eta:  2m | <1% bad
        warmup ===============                        1150/3000 | eta:  2m | <1% bad
        warmup ===============                        1200/3000 | eta:  2m | <1% bad
        warmup ================                       1250/3000 | eta:  2m | <1% bad
        warmup ================                       1300/3000 | eta:  2m | <1% bad
        warmup =================                      1350/3000 | eta:  2m | <1% bad
        warmup ==================                     1400/3000 | eta:  2m | <1% bad
        warmup ==================                     1450/3000 | eta:  2m | <1% bad
        warmup ===================                    1500/3000 | eta:  2m | <1% bad
        warmup ====================                   1550/3000 | eta:  1m | <1% bad
        warmup ====================                   1600/3000 | eta:  1m | <1% bad
        warmup =====================                  1650/3000 | eta:  1m | <1% bad
        warmup ======================                 1700/3000 | eta:  1m | <1% bad
        warmup ======================                 1750/3000 | eta:  1m | <1% bad
        warmup =======================                1800/3000 | eta:  1m | <1% bad
        warmup =======================                1850/3000 | eta:  1m | <1% bad
        warmup ========================               1900/3000 | eta:  1m | <1% bad
        warmup =========================              1950/3000 | eta:  1m | <1% bad
        warmup =========================              2000/3000 | eta:  1m | <1% bad
        warmup ==========================             2050/3000 | eta:  1m | <1% bad
        warmup ===========================            2100/3000 | eta:  1m | <1% bad
        warmup ===========================            2150/3000 | eta:  1m | <1% bad
        warmup ============================           2200/3000 | eta:  1m | <1% bad
        warmup ============================           2250/3000 | eta: 48s | <1% bad
        warmup =============================          2300/3000 | eta: 45s | <1% bad
        warmup ==============================         2350/3000 | eta: 42s | <1% bad
        warmup ==============================         2400/3000 | eta: 39s | <1% bad
        warmup ===============================        2450/3000 | eta: 36s | <1% bad
        warmup ================================       2500/3000 | eta: 33s | <1% bad
        warmup ================================       2550/3000 | eta: 29s | <1% bad
        warmup =================================      2600/3000 | eta: 26s | <1% bad
        warmup ==================================     2650/3000 | eta: 23s | <1% bad
        warmup ==================================     2700/3000 | eta: 20s | <1% bad
        warmup ===================================    2750/3000 | eta: 16s | <1% bad
        warmup ===================================    2800/3000 | eta: 13s | <1% bad
        warmup ====================================   2850/3000 | eta: 10s | <1% bad
        warmup =====================================  2900/3000 | eta:  7s | <1% bad
        warmup =====================================  2950/3000 | eta:  3s | <1% bad
        warmup ====================================== 3000/3000 | eta:  0s | <1% bad
    ## 
      sampling                                           0/1000 | eta:  ?s          
      sampling ==                                       50/1000 | eta: 38s          
      sampling ====                                    100/1000 | eta: 39s          
      sampling ======                                  150/1000 | eta: 45s          
      sampling ========                                200/1000 | eta: 40s          
      sampling ==========                              250/1000 | eta: 38s          
      sampling ===========                             300/1000 | eta: 37s          
      sampling =============                           350/1000 | eta: 35s          
      sampling ===============                         400/1000 | eta: 33s          
      sampling =================                       450/1000 | eta: 30s          
      sampling ===================                     500/1000 | eta: 28s          
      sampling =====================                   550/1000 | eta: 26s          
      sampling =======================                 600/1000 | eta: 22s          
      sampling =========================               650/1000 | eta: 20s          
      sampling ===========================             700/1000 | eta: 17s          
      sampling ============================            750/1000 | eta: 15s          
      sampling ==============================          800/1000 | eta: 12s          
      sampling ================================        850/1000 | eta:  9s          
      sampling ==================================      900/1000 | eta:  6s          
      sampling ====================================    950/1000 | eta:  3s          
      sampling ====================================== 1000/1000 | eta:  0s

``` r
summary(draws)
```

    ## 
    ## Iterations = 1:1000
    ## Thinning interval = 1 
    ## Number of chains = 4 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##          Mean        SD  Naive SE Time-series SE
    ## K     2.02046 0.0573854 9.073e-04      1.835e-03
    ## a     0.02569 0.0042425 6.708e-05      3.164e-04
    ## H     0.38097 0.0120478 1.905e-04      4.344e-04
    ## sigma 0.01990 0.0001989 3.144e-06      1.089e-05
    ## r     0.05489 0.0085158 1.346e-04      6.163e-04
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%     25%     50%     75%   97.5%
    ## K     1.91220 1.98132 2.01761 2.05675 2.14039
    ## a     0.01699 0.02291 0.02576 0.02879 0.03319
    ## H     0.35908 0.37285 0.38080 0.38859 0.40509
    ## sigma 0.01953 0.01976 0.01990 0.02004 0.02026
    ## r     0.03777 0.04911 0.05503 0.06119 0.07044

``` r
#bayesplot::mcmc_trace(draws)
```

``` r
samples <-  
  map_dfr(draws, 
          function(x) data.frame(x, t = 1:dim(x)[1]), 
          .id = "chain") %>% 
  gather(variable, value, -t, -chain)

samples %>%  
  ggplot(aes(t,value, col=chain, group=chain)) + 
  geom_line() +
  facet_wrap(~variable, scales = "free") + 
  scale_color_viridis_d()
```

![](greta_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
