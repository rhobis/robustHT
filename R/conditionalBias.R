#' Compute conditional bias
#'
#' Compute bias of the Horvitz-Thompson estimator conditionally to sample units
#'
#' @param y vector of observations of the variable of interest
#' @param pk vector of first-order inclusion probabilities; used only with
#' \code{sampling="ups"} and \code{sampling="poisson"}
#' @param pkl matrix of joint-inclusion probabilities; used only with \code{sampling="ups"}
#' @param n integer, the sample size; used only with \code{sampling="srs"}
#' @param N integer, the population size; used only with \code{sampling="srs"}
#' @param sample boolean, should the sample estimator be computed? Ignored for \code{sampling="poisson"}
#' @param sampling string indicating whether the conditional bias is to be computed
#' for a generic unequal probability sampling design (\code{"ups"}), a simple
#' random sampling (\code{"srs"}) or for Poisson sampling (\code{"poisson"}).
#'
#' @return a numeric vector of conditional bias



cond_bias <- function(y, pk, pkl, n, N, sample=TRUE, sampling = c("ups", "srs", "poisson")){

    sampling <- match.arg(sampling)

    ### check input

    if( sampling=="ups" & ( missing(pk) | missing(pkl) ) ) stop("pk and pkl objects are required for UPS sampling")
    if( sampling=="srs" & ( missing(n) | missing(N) ) ) stop("n and N must be specified for Simple Random Sampling")
    if( sampling=="poisson" & missing(pk) ) stop("pk object is required for Poisson sampling")
    if( sampling %in% c("ups", "poisson") ){
        if( length(y) != length(pk) ) stop("y and pk must have the same length")
    }
    if( sampling=="srs" ){
        if( sample & length(y)==N ){
            message( paste0("You chose to compute sample estimate on population data, ",
                            "exact formula for population conditional bias will be used instead") )
            sample <- FALSE
        }
        if( !( length(y) %in% c(n, N) ) ){
            stop( paste0("y vector length should be equal to the size of either ",
                         "the sample or the population") )
        }
    }

    ### conditional bias
    if( sampling == "ups" ){
        if( sample ){
            pp   <- outer(pk, pk, '*')
            b    <- (pkl-pp) / (pk*pkl)  * y
            b    <- apply(b,2,sum) #sample estimate of conditional bias for HT estimator
        }else{
            pp   <- outer(pk, pk, '*')
            b    <- (pkl-pp) / pp  * y
            b    <- apply(b,2,sum) #conditional bias for HT estimator
        }
    }else if( sampling == "srs" ) {
        if( sample ){
            b <- ( n/(n-1) ) * (N/n - 1) * ( y-mean(y) )
        }else{
            b <- ( N/(N-1) ) * ( N/n - 1 ) * ( y - mean(y) )
        }
    }else if (sampling == "poisson" ){
        b <- ( 1/pk - 1 ) * y
    }

    return( b )
}

