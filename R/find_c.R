### ----------------------------------------------------------------------------

#' Find the minimum c value iteratively
#'
#' Compute the value c that minimizes the maximum absolute conditional Bias of
#' the robust HT estimator over a grid of values.
#'
#' @param cBias vector of conditional bias for a sample
#'
#' @return a scalar, representing the value c that minimizes the maximum Bias
#'     of the robust HT estimator over a grid of values.
#'
#' @export


find_cmin_iter <- function(cBias, niter=1000){

    if(length(cBias)<1) return( NULL )

    acb    <- abs(cBias)
    rb     <- range(acb)
    cval   <- seq(rb[1],rb[2], length=niter)
    deltaC <- lapply(cval,
                     function(cv){
                         ## Huber function: sign(z)*min(abs(z), c)
                         psiB <- ifelse( acb < cv, cBias, sign(cBias)*cv)
                         return( sum(psiB-cBias) )
                     } )
    #estimate of conditional bias for the robust HT estimator
    b_rht <- lapply(deltaC,
                    function(D){
                        cBias + D
                    } )
    maxB <- sapply(b_rht, function(B) max( abs(B) ))
    ind_min <- which(maxB == min(maxB, na.rm = TRUE))
    cmin <- cval[ ind_min ]
    cmin <- max( cmin )

    return( cmin )
}


