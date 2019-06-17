greta.Rmd
================
Carl Boettiger
5/30/2019

``` r
library(tidyverse)
library(greta) # remotes::install_github("greta-dev/greta")
tensorflow::use_session_with_seed(42)
#set.seed(123456)
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
wide <- data %>% spread(reps, x) %>% select(-t) %>% 
  select(2) %>% # which / how many reps we use
  as.matrix()
n <- dim(wide)[1]
```

``` r
x_t1 <- wide[-1,] # data$x[-1]
x_t <- wide[-n,] # data$x[-100]
dim(x_t)
```

    ## NULL

``` r
dim(x_t1)
```

    ## NULL

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
system.time({
  draws <- mcmc(m, n_samples = 1000, warmup = 3000, chains = 4, verbose = FALSE)
})
```

    ##    user  system elapsed 
    ##  92.295  17.871  38.647

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
    ## K     2.13812 0.1208677 1.911e-03      9.658e-03
    ## a     0.04405 0.0108523 1.716e-04      1.142e-03
    ## H     0.41827 0.0280719 4.439e-04      2.503e-03
    ## sigma 0.02016 0.0004642 7.340e-06      1.243e-05
    ## r     0.08834 0.0225332 3.563e-04      2.756e-03
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%     25%     50%     75%   97.5%
    ## K     1.96152 2.06334 2.12108 2.18861 2.46614
    ## a     0.01758 0.03777 0.04608 0.05269 0.05825
    ## H     0.37015 0.40158 0.41729 0.43324 0.48367
    ## sigma 0.01928 0.01986 0.02015 0.02047 0.02111
    ## r     0.03312 0.07543 0.09274 0.10555 0.12124

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

``` r
#Q = 5
true <- data.frame(a = 0.023, r = .05, K = 2, H = .38, sigma = .02) %>%
  gather(variable, value)
```

``` r
samples %>% ggplot() + 
  geom_histogram(aes(value), bins = 30)  +
  geom_vline(data = true, aes(xintercept = value), col = "red", lwd = 1) + 
  facet_wrap(~variable, scales = "free")
```

![](greta_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
