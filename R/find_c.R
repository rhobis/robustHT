### ----------------------------------------------------------------------------

#' Find the minimum c value over a grid of values
#'
#' Compute the value c that minimizes the maximum absolute conditional Bias of
#' the robust HT estimator over a grid of values.
#'
#' @param cBias vector of conditional bias for a sample
#' @param ngrid integer scalar indicating the length of the grid of values that
#'         is generated to estimate the optimum value for c
#'
#' @return A scalar, representing the value c that minimizes the maximum Bias
#'     of the robust HT estimator over a grid of values.
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
#' # Find the minimum c
#' find_cmin(cb, ngrid = 200)
#'
#'
#' @export
#'


find_cmin <- function(cBias, ngrid=4*length(cBias)){

    if(!is.numeric(cBias) | !is.vector(cBias) | is.list(cBias) )
        stop("Argument cBias must be a vector!")
    if(length(cBias)<2) stop("Argument cBias must have length >2")


    if( is.array(ngrid) | is.list(ngrid) | is.data.frame(ngrid) | !is.numeric(ngrid) )
        stop("Argument ngrid must be numeric scalar!")
    if( length(ngrid)>1) ngrid <- ngrid[1]
    ngrid <- ceiling(ngrid)

    acb    <- abs(cBias)
    rb     <- range(acb)
    cval   <- seq(rb[1],rb[2], length=ngrid)
    deltaC <- lapply(cval,
                     function(cv){
                         ## Huber function: sign(z)*min(abs(z), c)
                         psi <- sign(cBias) * pmin(acb, cv)
                         return( sum(psi-cBias) )
                     } )
    #conditional bias estimates for the robust HT estimator
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


