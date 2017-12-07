
## LagBasis
## -------------------------------------------------------------------
#' @title Create and Manipulate Lag Basis Functions
#' @inherit basis references
#'
#' @description
#' S4 class object to store and query components of
#' lag basis functions. User interface for creating this class can be found in
#' \code{\link{basis}}.
#'
#' @slot x
#'   original lag data
#' @slot x.center
#'   store the value the lag data was centered to. In theory
#'   this is useful for the non-yet-implemented \code{predict} method
#' @slot x.scale
#'   store the value the lag data was scaled by. Again, should be
#'   useful primarily for the \code{predict} method in the future
#' @slot C0
#'   C_0 part of basis matrix
#' @slot K1
#'   K_1 part of basis matrix
#' @slot dist.fun
#'   whatever distance function was used to compute the pointwise
#'   distances between the elements of \code{x}
#'
#' @details
#' See \code{\link{basis}} for details of the decomposition
#'
#' @name LagBasis
#'
LagBasis <- setClass("LagBasis",
         slots = c(
           x  = "numeric",
           x.center = "numeric",
           x.scale = "numeric",
           C0 = "matrix",
           K1 = "matrix",
           dist.fun = "function"
         ),
         prototype = list(
           x = NA_real_,
           x.center = 0,
           x.scale = 1,
           C0 = matrix(),
           K1 = matrix()
         ))




## SmoothLag
## -------------------------------------------------------------------
#' @title Lag Matrix With Applied Smoothing
#'
#' @description
#' An S4 class object for representing a smoothed lag matrix in addition
#' to containing details about the basis smoothing. Intended for use within
#' \code{\link{dlm}} modeling framework to assist extraction of basis
#' components treated as "fixed" and "random" effects. Inherits from
#' \code{matrix}.
#'
#' @slot basis
#'   A \code{\link{LagBasis}} smoothing object containing details about the
#'   lag and the smoothing parameters used
#'
#' @slot .Data
#'   Contains the "fixed effects" components of the smoothed lag function.
#'   This scheme is intended to work conveniently with
#'   \code{stats::model.matrix}
#'
#' @slot random
#'   Contains the "random effects" components of the smoothed lag function
#'
#' @slot terms
#'   Character vector contianing the function name and deparsed arguments
#'   from whatever smoothing function generated the \code{SmoothLag} object
#'
#' @name SmoothLag
#'
SmoothLag <- setClass("SmoothLag",
                       slots = c(
                         basis = "LagBasis",
                         random = "matrix",
                         signature = "character"
                       ),
                       contains = "matrix"
                       )




## dlMod
## -------------------------------------------------------------------
#' @title Distributed lag models
#'
#' @slot resp An \code{\link[lme4]{lmResp}} object to store a
#' mixed-model response variable
#'
#' @slot bases A list of \code{\link{LagBasis}} objects corresponding
#' to the unique set of lag bases used to fit the model
#'
#' @slot index A (named) integer vector providing the index of the
#' basis set in \code{bases} corresponding to each lag term in the
#' model
#'
#' @name dlMod
#'
dlMod <- setClass("dlMod",
                  slots = c(
                    resp = "lmResp",   # from lme4
                    bases = "list",
                    index = "numeric"
                  ),
                  contains = "merMod"  # also from lme4
                  )

## ldlMod <- setClass("ldlMod", contains = "dlMod")
## gdlMod <- setClass("gdlMod", contains = "dlMod")



