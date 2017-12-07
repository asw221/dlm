
## is.LagBasis <- function(object, ...)  inherits(object, "LagBasis")
## is.SmoothLag <- function(object, ...) inherits(object, "SmoothLag")
## is.dlMod <- function(object, ...)     inherits(object, "dlMod")


makeDlMod <- function(object, ...) UseMethod("makeDlMod", object)



## cholfVar
## -------------------------------------------------------------------
#' @title Extract Cholesky factor of inverse Information matrix
#' @name cholfVar
setGeneric("cholfVar", function(object, ...) standardGeneric("cholfVar"))


## Sigma
## -------------------------------------------------------------------
#' @title Extract model coefficients variance matrix
#' @name Sigma
setGeneric("Sigma", function(object, ...) standardGeneric("Sigma"))


## vcoef
## -------------------------------------------------------------------
#' @title Vectorized coefficients
#'
#' @description
#' Extract fixed and random effects coefficient vector,
#' \eqn{(\beta, b)'} from a fitted \code{\link{dlMod}} object
#'
#'
#' @return A numeric vector: with variable names in the case of \code{vcoef}
#'
#' @name vcoef
setGeneric("vcoef", function(object, ...) standardGeneric("vcoef"))
setGeneric("vcoef0", function(object, ...) standardGeneric("vcoef0"))


## lagIndex
## -------------------------------------------------------------------
#' @title Extract list of indices of lag terms
#' @name lagIndex
setGeneric("lagIndex", function(object, ...) standardGeneric("lagIndex"))



## scaleMat
## -------------------------------------------------------------------
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
#' @name scaleMat
setGeneric("scaleMat",
           function(object, ...) standardGeneric("scaleMat")
           )




## changePoint
## -------------------------------------------------------------------
#' @title Lag Coefficient Change Points
#'
#' @return An integer vector
#' @name changePoint
setGeneric("changePoint", function(object, ...) standardGeneric("changePoint"))







## omega
## -------------------------------------------------------------------
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
#'
#' @name omega
setGeneric("omega", function(object, ...) standardGeneric("omega"))


