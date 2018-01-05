

## omega
## -------------------------------------------------------------------
#' @describeIn omega Method for \code{"\link{LagBasis}"} objects
setMethod("omega", signature = "LagBasis",
          function(object, ...)  cbind(object@C0, object@K1)
          )



## predict.LagBasis
## -------------------------------------------------------------------
#' @title Predict new values for fitted lag basis
#'
#' @param object A \code{\link{LagBasis}} object
#'
#' @description
#' Will update w/ description. Not yet implemented
#'
predict.LagBasis <- function(object, x, ...) {
  stop ("Not yet Implemented")
  cx <- (object@x - object@x.center) / object@x.scale
  dcmp <- .ss.decomp(cx, object@dist.fun, ...)
  cx.new <- (x - object@x.center) / object@x.scale
  C0 <- cbind(1, cx.new)

}
## predict.LagBasis

