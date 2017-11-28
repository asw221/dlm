

basis <- function(x, center = TRUE, scale = FALSE, ...,
                  .fun = function(x, ...) abs(outer(x, x, "-"))^3
                  ) {
    ## scale and center lag, if desired
  cntr <- 0
  scl <- 1
  if (center && is.logical(center))
    cntr <- mean(x)
  else if (center)
    cntr <- center

  if (scale && is.logical(scale))
    scl <- sd(x)
  else if (scale)
    scl <- scale
  cx <- c(x - cntr) / scl

  ## compute basis functions
  C0 <- cbind(1, cx)
  colnames (C0) <- x[1:NCOL(C0)]
  C1 <- .fun(cx, ...)
  M1 <- qr.Q(qr(cbind(C0, C1)))[, -(1:2)]
  S <- svd(t(M1) %*% C1 %*% M1)
  M2.inv.sqrt <- S$v %*% diag(1 / sqrt(S$d)) %*% t(S$u)
  K1 <- C1 %*% M1 %*% M2.inv.sqrt
  ## rescl <- max(max(abs(K1)), max(abs(C0)))
  colnames (K1) <- x[-(1:NCOL(C0))]
  LagBasis(x = x, x.center = cntr, x.scale = scl,
           C0 = C0, K1 = K1
           ## C0 = C0 / rescl, K1 = K1 / rescl
           )
}





#' @title Cubic Radial Basis
#'
#' @description
#' Construct a natural cubic radial basis function for a given
#' distance set vector and apply as a linear transformation of a
#' covariate matrix.
#'
#' @param x
#'   a vector of values to construct the basis from. Missing values
#'   are not allowed.
#'
#' @param Z
#'   a covariate matrix (or object that can be coerced to a \code{matrix})
#'   to take the basis transformation of.
#'   \code{length(x)} should be the same as \code{ncol(Z)}.
#'   Missing values are not allowed.
#'
#' @param ...
#'   arguments to be passed to \code{\link{lag.basis}}
#'
#' @details
#'   At the time of writing, \code{cr} is little more than a convenient
#'   wrapper to \code{\link{lag.basis}} intended to simplify the task of
#'   specifying lag terms in a model \code{formula}. The longer term goal
#'   is for \code{cr} (and potentially other basis functions) to replace
#'   \code{lag.basis} entirely.
#'
#' @return
#' An S4 object of class \code{\link{SmoothBasis}}.
#'
## basis function extensions should be of class "SmoothLag"

cr <- function(x, Z, ...) {
  if (any(is.na(x)) || any(is.na(Z)))
    stop ("missing values not allowed")
  if (length(x) != ncol(Z <- as.matrix(Z)))
    stop ("arguments do not have compatible dimensions")
  B <- basis(x, ...)
  SmoothLag(Z %*% B@C0, random = Z %*% B@K1,
            basis = B,
            signature = deparse(sys.call())
            )
}


