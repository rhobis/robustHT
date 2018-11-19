robustHT
======================================================

[![Travis-CI Build Status](https://travis-ci.org/rhobis/UPSvarApprox.svg?branch=master)](https://travis-ci.org/rhobis/UPSvarApprox)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/robustHT)](https://cran.r-project.org/package=robustHT)
[![](https://cranlogs.r-pkg.org/badges/grand-total/robustHT)](https://cran.r-project.org/package=robustHT)


Description 
-----------------

UPSvarApprox provides functions for the approximation of the variance of the 
Horvitz-Thompson total estimator in Unequal Probability Sampling
using only first-order inclusion probabilities.

The main functions are:

- `cond_bias()`: computes or estimates the conditional bias; 
- `RHTest()`: computes the Robust Horvitz-Thompson estimator;



Installation
------------

The development version of the package can be installed from GitHub:

``` r
# if not present, install 'devtools' package
install.packages("devtools")
devtools::install_github("rhobis/UPSvarApprox")
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
- For more information, please contact the manteiner at `roberto.sichera@unipa.it`. 
