

## omega
## -------------------------------------------------------------------
#' @describeIn omega Method for "\code{LagBasis}" objects
setMethod("omega", signature = "LagBasis",
          function(object, ...)  cbind(object@C0, object@K1)
          ## function(object, ...)  cbind(object@K1, object@C0)
          )



## predict.LagBasis
## -------------------------------------------------------------------
#' @title Predict New Values for Fitted Lag Basis
#'
#' @param object A \code{\link{LagBasis}} object
#'
#' @description
#' Not yet implemented
#'
predict.LagBasis <- function(object, ...) {
  stop ("Not yet Implemented")
}
## predict.LagBasis

