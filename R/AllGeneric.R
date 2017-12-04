## New Generics
## -------------------------------------------------------------------


is.LagBasis <- function(object, ...)  inherits(object, "LagBasis")
is.SmoothLag <- function(object, ...) inherits(object, "SmoothLag")
is.dlMod <- function(object, ...)     inherits(object, "dlMod")


makeDlMod <- function(object, ...) UseMethod("makeDlMod", object)




## random effects indices (list)
setGeneric("reIndex", function(object, ...) standardGeneric("reIndex"))

## Variance matrix for regression coefficients
setGeneric("Sigma", function(object, ...) standardGeneric("Sigma"))


## Extract Cholesky factor of inverse Information matrix
setGeneric("cholfVar", function(object, ...) standardGeneric("cholfVar"))

## Extract [ranef, fixef] coefficient vector
setGeneric("vcoef", function(object, ...) standardGeneric("vcoef"))
setGeneric("vcoef0", function(object, ...) standardGeneric("vcoef0"))


## Extract list of indices of lag terms
setGeneric("lagIndex", function(object, ...) standardGeneric("lagIndex"))



#' @title Extract Lag Basis Matrix
#'
#' @description
#' Extract lag basis matrix, \eqn{\Omega}
#'
#' @return
#' A square numeric matrix
#'
#' @details
#' \code{omega} is S4 generic.
setGeneric("omega", function(object, ...) standardGeneric("omega"))






#' @title Extract Standard Errors of Fixed Effects
#'
#' @description
#' Returns a vector of fixed effects standard errors from a fitted
#' mixed effects model.
#'
#' @return
#' A numeric vector
#'
#' @details
#' \code{fixef} is S4 generic.
#'
setGeneric("se.fixef",
           function(object) standardGeneric("se.fixef")
           )


#' @title Extract Standard Errors of Random Effects
#'
#' @description
#' Returns a vector of random effects standard errors from a fitted
#' mixed effects model.
#'
#' @return
#' A numeric vector
#'
#' @details
#' \code{fixef} is S4 generic.
#'
setGeneric("se.ranef",
           function(object) standardGeneric("se.ranef")
           )



#' @title Extract Distributed Lag Scale Matrix
#'
#' @description
#' Return lag coefficient scale matrix, S, such that the distributed lag
#' coefficients fit by the model are obtained via the transformation
#' \eqn{\beta = S \theta}{\beta = S * \theta}. S should be invertable.
#'
#' @param object a fitted model object
#'
#' @return A square numeric matrix
#'
#' @details
#' \code{scaleMat} is S4 generic.
#'
setGeneric("scaleMat",
           function(object, ...) standardGeneric("scaleMat")
           )





#' @title Lag Coefficient Change Points
#'
#' @return An integer vector
#'
setGeneric("changePoint", function(object, ...) standardGeneric("changePoint"))












## Methods for LagBasis class
## -------------------------------------------------------------------


#' @describeIn omega Method for "\code{LagBasis}" objects
setMethod("omega", signature = "LagBasis",
          function(object, ...)  cbind(object@C0, object@K1)
          ## function(object, ...)  cbind(object@K1, object@C0)
          )


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


## Methods for LagBasis class
## -------------------------------------------------------------------


#' @describeIn omega Method for "\code{LagBasis}" objects
setMethod("omega", signature = "LagBasis",
          function(object, ...)  cbind(object@C0, object@K1)
          )


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
