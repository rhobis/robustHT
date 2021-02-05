#' Compute Delta(cmin)
#'
#' Compute Delta(cmin) quantity for the Robust Horvitz-Thompson estimator
#'
#'
#' @param cBias vector of conditional bias, computed by function \code{conditional_bias()}
#'
#' @return a scalar, representing the Delta value corresponding to the minimum c
#'
#' @examples
#' # Generate population data
#' N <- 50; n <- 5
#'
#' set.seed(0)
#' x <- rgamma(500, scale=10, shape=5)
#' y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#'
#' # Select sample
#' pik <- n * x/sum(x)
#' s   <- sample(N, n)
#' ys <- y[s]
#' piks <- pik[s]
#'
#' # Compute conditional bias
#' cb <- conditional_bias(y=ys, pk=piks, sampling = "poisson")
#'
#' # Find Delta(cmin)
#' Delta_cmin(cb)
#'
#' @export

Delta_cmin <- function(cBias){
    rangeB <- range(cBias)
    D <- -0.5 * (rangeB[1] + rangeB[2])
    return(D)
}
