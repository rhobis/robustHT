### ----------------------------------------------------------------------------

#' Robust Horvitz-Thompson estimator
#'
#' Estimate the Robust Horvitz-Thompson by
#' Beaumont, Haziza, Ruiz-Gazen (2013)
#'
#' @param ys vector of sample observations
#' @param method string that determines how to compute the penalised part of
#'        the robust estimator. Available methods are 'find_c' and 'Delta_min',
#'        see section Details for more information.
#' @param pks vector of sample first-order inclusion probabilities; only required if
#'         \code{sampling="ups"} or \code{sampling="poisson"}
#' @param pkls matrix of size nxn with joint-inclusion probabilities;
#'         required only with \code{sampling="ups"}
#' @param n integer, the sample size; used only with \code{sampling="srs"}
#' @param N integer, the population size; used only with \code{sampling="srs"}
#' @param sampling string indicating whether the conditional bias is to be computed
#' for a generic unequal probability sampling design (\code{"ups"}), a simple
#' random sampling (\code{"srs"}) or for Poisson sampling (\code{"poisson"}).
#' @param grid_length  integer scalar indicating the length of the grid of values that
#'         is generated to estimate the optimum value for c
#'
#' @details
#' The argument \code{method} accepts as input the strings 'find_c' or 'Delta_min'.
#' If the former is chosen, the \eqn{Delta}{\Delta(c)} component in the robust estimator
#' is computed by finding the optimal \code{c} value over a grid of values, while
#' the former computes the \eqn{Delta}{\Delta(c)} corresponding to the optimum \code{c}
#' by using the formula provided by Beaumont, Haziza and Ruiz-Gazen (2013), section 3.3
#'
#'
#'
#' @return The Robust Horvitz-Thompson estimate
#'
#'
#' @examples
#' ### Generate population data ---
#' N <- 50; n <- 5
#'
#' set.seed(0)
#' x <- rgamma(500, scale=10, shape=5)
#' y <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#'
#' pik <- n * x/sum(x)
#' s   <- sample(N, n)
#' ys <- y[s]
#' piks <- pik[s]
#'
#' RHTestimator(ys, piks, method='find_c', sampling='poisson')
#' RHTestimator(ys, piks, method='Delta_min', sampling='poisson')
#'
#'
#'
#' @export
#'


RHTestimator <- function(ys, pks=NULL, method=c('find_c', 'Delta_min'),
                        sampling = c("ups", "srs", "poisson"),
                        pkls=NULL, n=NULL, N=NULL,
                        grid_length = 4*length(ys) ){


    ## Check input
    method <- match.arg(method, c('find_c', 'Delta_min'))
    sampling <- match.arg(sampling, c("ups", "srs", "poisson"))

    if(!is.numeric(ys) | !is.vector(ys) | is.list(ys))
        stop("Argument ys must be a numeric vector!")
    if(length(ys)<2) stop("Argument ys must have length > 1")
    if( any(ys %in% c(NA, NaN)) ) stop("Vector ys includes NA or NaN values!")

    if(method=='find_c' & (grid_length %in% c(NaN, NA) | is.null(grid_length)))
        stop('Invalid value for argument grid_length!')

    if(sampling == 'srs' & (missing(n) | missing(N)) ){
        stop("Both arguments n and N are required when sampling='srs' ")
    }else{
        if(!is.numeric(pks) | !is.vector(pks) | is.list(pks))
            stop("Argument pks must be a numeric vector!")
        if(length(pks)<2) stop("Argument pks must have length > 1")
        if( any(pks %in% c(NA, NaN)) ) stop("Vector pks includes NA or NaN values!")

        if( length(ys) != length(pks) ) stop("Arguments ys and pks must have the same length!")
    }

    if(sampling == 'ups'){
        if( missing(pkls) ) stop("Argument pkls is required with sampling='ups' ")
        if(!is.matrix(pkls)) stop("Argument pkls must be a square matrix!")
        if( diff(dim(pkls)) != 0 ) stop("Argument pkls must be a square matrix!")
        if( any(is.na(pkls)) ) stop("The pkls matrix includes NA values! ")
        if( any(is.nan(pkls)) ) stop("The pkls matrix includes NaN values! ")
        if( nrow(pkls) != length(ys)) stop(paste0("The number of rows and columns ",
                                                  "of argument pkl must be equal to the length of ys") )
    }

    ## Estimate conditional bias
    CB <-
        if(sampling=='srs'){
            conditional_bias(y=ys, n=n, N=N, sample=TRUE, sampling = "srs")
        }else if(sampling == 'poisson'){
            conditional_bias(y=ys, pk=pks, sampling = "poisson")
        }else{
            conditional_bias(y=ys, pk=pks, pkl=pkls, sample=TRUE, sampling = "ups")
        }


    ## Compute Delta
    if(method == 'Delta_min'){
        Delta <- Delta_cmin(CB)
    }else{
        cmin  <- find_cmin(CB, ngrid = grid_length)
        psi   <- sign(CB) * pmin(abs(CB), cmin)
        Delta <- sum( psi - CB )
    }

    ## Compute the RHT estimator
    HTtot <-
        if(sampling == 'srs'){
            N*mean(ys)
        }else sum(ys/pks)

    return( HTtot + Delta )


}


### ----------------------------------------------------------------------------

#' Compute Delta(cmin) quantity for the Robust Horvitz-Thompson estimator
#'
#'
#' @param cBias vector of conditional bias, computed by function \code{conditional_bias()}
#'
#' @return a scalar, representing the Delta value corresponding to the minimum c
#'
#' @export

Delta_cmin <- function(cBias){
    rangeB <- range(cBias)
    D <- -0.5 * (rangeB[1] + rangeB[2])
    return(D)
}




