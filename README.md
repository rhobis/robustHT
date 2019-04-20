robustHT
======================================================

[![Travis-CI Build Status](https://travis-ci.org/rhobis/robustHT.svg?branch=master)](https://travis-ci.org/rhobis/robustHT)


Description 
-----------------

This is a simple package for obtaining the conditional bias of the Horvitz-Thompson
estimator and for estimating the Robust Horvitz-Thompson Estimator (Beaumont, Haziza and Ruiz-Gazen, 2013).

The main functions are:

- `conditional_bias()`: computes or estimates the conditional bias; 
- `RHTestimator()`: estimates the Robust Horvitz-Thompson total;



Installation
------------

The development version of the package can be installed from GitHub:

``` r
# if not present, install 'devtools' package
install.packages("devtools")
devtools::install_github("rhobis/robustHT")
```

Usage
-----

``` r
library(robustHT)

### Generate population data ---
N <- 50; n <- 5

set.seed(0)
x <- rgamma(500, scale=10, shape=5)
y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )

pik <- n * x/sum(x)
s   <- sample(N, n)

ys <- y[s]
piks <- pik[s]


### Estimate conditional bias
cb <- conditional_bias(ys, piks, sampling='poisson')


### Estimate the Robust Horvitz-Thompson total
RHTestimator(ys, piks, method='find_c', sampling='poisson', grid_length=10000)
RHTestimator(ys, piks, method='Delta_min', sampling='poisson')

```






More
----

- Please, report any bug or issue [here](https://github.com/rhobis/robustHT/issues).
- For more information, please contact the maintainer at `roberto.sichera@unipa.it`. 
