#' Compute the conditional bias
#'
#' Compute the bias of the Horvitz-Thompson estimator conditionally to sample units
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
#'
#' @export



conditional_bias <- function(y, pk=NULL, pkl=NULL, n=NULL, N=NULL, sample=TRUE,
                             sampling = c("ups", "srs", "poisson")){

    sampling <- match.arg(sampling, c("ups", "srs", "poisson"))

    ### Check input

    if(!is.numeric(y) | !is.vector(y) | is.list(y))
        stop("Argument y must be a numeric vector!")
    if(length(y)<2) stop("Argument y must have length > 1")
    if( any(y %in% c(NA, NaN)) ) stop("Vector y includes NA or NaN values!")

    if( sampling=="ups" & ( missing(pk) | missing(pkl) ) )
        stop("Arguments pk and pkls are required for sampling='ups' ")
    if( sampling=="srs" & ( missing(n) | missing(N) ) )
        stop("Arguments n and N must be specified for sampling='srs' ")
    if( sampling=="poisson" & missing(pk) )
        stop("Argument pk is required for sampling='poisson' ")

    if( sampling=="srs" ){
        if( missing(n) | missing(N) )
            stop("Both arguments n and N are required when sampling='srs' ")

        if( sample & length(y)==N ){
            message( paste0("You chose to compute sample estimate on population data, ",
                            "exact formula for population conditional bias will be used instead") )
            sample <- FALSE
        }
        if( !( length(y) %in% c(n, N)) )
            stop( paste0("y vector length should be equal to the size of either ",
                         "the sample or the population") )
    } else {
        if(!is.numeric(pk) | !is.vector(pk) | is.list(pk))
            stop("Argument pk must be a numeric vector!")
        if(length(pk)<2) stop("Argument pk must have length > 1")
        if( any(pk %in% c(NA, NaN)) ) stop("Vector pk includes NA or NaN values!")
        if( length(y) != length(pk) ) stop("Arguments y and pk must have the same length!")
    }

    if(sampling == 'ups'){
        if( missing(pkl) ) stop("Argument pkl is required with sampling='ups' ")
        if(!is.matrix(pkl)) stop("Argument pkl must be a square matrix!")
        if( diff(dim(pkl)) != 0 ) stop("Argument pkl must be a square matrix!")
        if( any(is.na(pkl)) ) stop("The pkl matrix includes NA values! ")
        if( any(is.nan(pkl)) ) stop("The pkl matrix includes NaN values! ")
        if( nrow(pkl) != length(y)) stop(paste0("The number of rows and columns ",
                                                "of argument pkl must be equal to the length of y") )
    }



    ### Compute conditional bias
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
            b <- ( n/(n-1) ) * (N/n - 1) * ( y - mean(y) )
        }else{
            b <- ( N/(N-1) ) * ( N/n - 1 ) * ( y - mean(y) )
        }
    }else if (sampling == "poisson" ){
        b <- ( 1/pk - 1 ) * y
    }

    return( b )
}

