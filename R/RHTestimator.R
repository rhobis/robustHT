

### ----------------------------------------------------------------------------

#' Robust Horvitz-Thompson estimator
#'
#' Estimate the Robust Horvitz-Thompson estimator according to
#' Beaumont, Haziza, Ruiz-Gazen (2013)
#'
#' @param ys vector of sample observed values
#' @param pks vector of sample inclusion probabilities
#' @param cBias vector of conditional bias, used if \code{Delta_cmin} is missing
#' @param c tuning parameter to compute robust estimate, used if Delta_cmin is missing
#' @param Delta_cmin Delta value, corresponding to the minimum value of the tuning
#' constant c, computed by function \code{min_Delta()}. Either this or the combination
#' of \code{cBiad} and \code{c} should be used.
#'
#' @return The Robust Horvitz-Thompson estimate
#'

RHTestimate <- function(ys, pks, cBias, cmin, Delta_cmin){

    cond1 <- !missing(Delta_cmin)
    cond2 <- ( !missing(cBias) & !missing(cmin))

    if( cond1 ){
        if( !is.numeric(Delta_cmin) ) stop('Delta_cmin must be numeric!')
        rht <- sum(ys/pks) + Delta_cmin

    }else if( cond2 ){
        if( !is.numeric(cBias) | !is.numeric(cmin) ) stop(' cBias and cmin must be numeric!')
        psi <- pmax( -cmin, pmin(cmin, cBias) )
        rht <- sum(ys/pks) + ( sum(psi) - sum(cBias) )

    }else stop("Either Delta_cmin or both cBias and c must be provided, ",
               "and they must be numeric.")

    return( rht )


}


### ----------------------------------------------------------------------------

#' Compute Delta(cmin) for the Robust Horvitz-Thompson estimator
#'
#'
#' @param cBias vector of conditional bias, computed by function \code{cond_bias()}
#'
#' @return a scalar, representing the Delta value corresponding to the minimum c
#'

Delta_cmin <- function(cBias){
    rangeB <- range(cBias)
    D <- -0.5 * (rangeB[1] + rangeB[2])
    return(D)
}
