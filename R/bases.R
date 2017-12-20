
## basis function extensions should be of class "SmoothLag"

## basis
## -------------------------------------------------------------------
#' @title Basis vector sets
#'
#' @description
#' Construct a set of basis vectors based on the distances between
#' input points.
#'
#' Users should not typically interact with \code{basis} directly.
#' Typical usage relies on calling basis application functions, like
#' \code{\link{cr}} (e.g. in \code{\link{dlm}} model
#' formulas).
#'
#' @usage
#' basis(x, center = TRUE, scale = FALSE, .fun = NULL, ...)
#'
#' @param x
#'   a set of points to measure distances between; the resultant
#'   distance matrix will be decomposed into a set of basis vectors
#' @param center
#'   if \code{TRUE} (the default), parameter \code{x} will be mean
#'   centered prior to computing distances. Otherwise, if given a
#'   \code{numeric} value, \code{x} will be centered at \code{center}
#' @param scale
#'   if \code{TRUE} (default = \code{FALSE}), parameter \code{x}
#'   will be scaled by \code{sd(x)}. Otherwise, if given a
#'   \code{numeric} value, \code{x} will be scaled by \code{scale}
#' @param .fun
#'   a function to compute distances between the values in \code{x}.
#'   The default is to compute pairwise cubed absolute distances. See
#'   Details
#' @param ...
#'   other parameters passed to \code{.fun}
#'
#' @details
#' Alternative distance functions, \code{.fun}, may be specified, and
#' error checking on the user's choice of \code{.fun} is deliberately
#' missing. Proper candidates for \code{.fun} should return a
#' (\code{length(x)} by \code{length(x)}) matrix, however, and values
#' are typically non-negative.
#'
#' In addition, new distance function definitions should follow the idiom:
#' \code{
#'   function(x, y, ...) {
#'     if (missing(y))  y <- x
#'     ...
#'   }
#' }
#'
#' The default value of \code{.fun} computes cubic radial distance,
#' which amounts to \code{abs(outer(x, y, "-"))^3}; distance matrix
#' decomposition follows Rupert, Wand, and Carroll
#' (2003). In particular, once \code{x} and \code{.fun} are chosen,
#' define distance matrix \code{C_1 = .fun(x, ...)}, and let
#'
#' \deqn{C_0 = [1^{(n \times 1)}, x]}{C_0 = [1, x]}
#' \deqn{C_1 = Q R}{C_1 = Q * R}
#' \deqn{M_1 = Q_{-(1:2)}}{M_1 = Q[-(1:2)]}
#' \deqn{K_1 = C_1 M_1 (M_1^T C_1 M_1)^{-\frac{1}{2}}}{K_1 = C_1 * M_1 * (M_1' * C_1 * M_1)^-0.5}
#'
#' where \eqn{A_{-j}}{A[-j]} denotes a matrix \eqn{A} with column(s) \eqn{j}
#' removed. Then distributed lag effects of interest are scaled by the
#' matrix \eqn{\Omega = [C_0, K_1]}{\Omega = [C_0, K_1]}
#'
#' @return An object of class \code{\link{LagBasis}}
#'
#' @seealso \code{\link{cr}}, \code{\link{dlm}}
#'
#' @references Rupert D, Wand MP, & Carroll RJ (2003) Semiparametric
#' Regression. New York: Cambridge University Press.
#'
#' @name basis
basis <- function(x, center = TRUE, scale = FALSE, .fun = NULL, ...) {
  ## setup .fun
  ## default to pairwise cubic absolute distance
  if (is.null(.fun))  .fun <- .cr
  else if (!is.function(.fun))
    stop (".fun must be a function")

  ## scale and center lag, if desired
  cntr <- 0
  scl <- 1
  if (center && is.logical(center))
    cntr <- mean(x, na.rm = TRUE)
  else if (center)
    cntr <- center

  if (scale && is.logical(scale))
    scl <- sd(x, na.rm = TRUE)
  else if (scale)
    scl <- scale
  cx <- c(x - cntr) / scl

  ## compute basis vectors
  ## if numerical issues come to the fore at some point,
  ## may consider norming the C0, K1 matrices
  dcmp <- .ss.decomp(cx, .fun, ...)
  colnames (dcmp$C0) <- x[1:NCOL(dcmp$C0)]
  colnames (dcmp$K1) <- x[-(1:NCOL(dcmp$C0))]
  LagBasis(x = x, x.center = cntr, x.scale = scl,
           C0 = dcmp$C0, K1 = dcmp$K1, dist.fun = .fun
           )
}
## basis





## Smoothing Help
## -------------------------------------------------------------------
#' @title Basis smoothing
#' @inherit basis references
#'
#' @description
#' Construct a set of basis vectors for distances between distributed
#' lag points, and apply as a linear transformation of a
#' concentration matrix.
#'
#' @param x
#'   a vector of values to construct the basis from. Missing values
#'   are not allowed
#' @param Z
#'   a covariate matrix (or object that can be coerced to a \code{matrix})
#'   to apply the linear basis transformation to.
#'   \code{length(x)} should be the same as \code{ncol(Z)}
#' @param ...
#'   arguments to be passed to \code{\link{basis}}
#'
#' @details
#'   These functions are little more than convenient
#'   wrappers to the function \code{\link{basis}} and the
#'   \code{\link{SmoothLag}} class constructor. They are intended to
#'   simplify the task of specifying lag terms in a model \code{formula}.
#'   The functions computes a set of basis vectors for
#'   parameter \code{x} and applies this basis as a linear transformation
#'   of the covariate/concentration matrix parameter, \code{Z}. For example,
#'   if \code{Z} is the identity matrix, the model fit will simply be
#'   the natural cubic spline of \code{x}.
#'
#'   Note that other basis extensions should always return an object
#'   that inherits from \code{\link{SmoothBasis}}
#'
#' @seealso \code{\link{basis}}, \code{\link{SmoothBasis}}
#'
#' @return
#' An S4 object of class \code{\link{SmoothBasis}}.
#'
#' @examples
#' ## load simulated data set and extract concentration matrix
#' data (simdata)
#' Conc <- simdata[, -(1:3)]  # First columns are Y, Age, and Gender
#'
#' ## radial lag (distance) each concentration was measured at
#' x <- seq(0.1, 10, length.out = ncol(Conc))
#' crb <- cr(x, Conc)
#'
#' @name smoothing
NULL




#' @describeIn smoothing natural cubic radial basis spline
cr <- function(x, Z, ...) {
  .Ignored(...)
  val <- sm(x, Z)
  val@signature <- deparse(sys.call())
  val
}
## cr


#' @describeIn smoothing user-defined smoothing
sm <- function(x, Z, ..., .fun = NULL) {
  if (any(is.na(x)) || any(is.na(Z)))
    stop ("missing values not allowed")
  if (is.data.frame(Z))  Z <- as.matrix(Z)
  if (length(x) != NCOL(Z))
    stop ("arguments do not have compatible dimensions")
  B <- basis(x, .fun = .fun, ...)
  SmoothLag(as.matrix(Z %*% B@C0), random = as.matrix(Z %*% B@K1),
            basis = B,
            signature = deparse(sys.call())
            )
}
## sm





.cr <- function(x, y, ...) {
  if (missing(y)) y <- x
  abs(outer(x, y, "-"))^3
}








## .ss.decomp
## -------------------------------------------------------------------
.ss.decomp <- function(x, .fun, ...) {
  ## decompose the basis vector set to control degrees of
  ## freedom in the downstream model. This bit follows
  ## Rupert Wand and Carroll (see references above)
  C0 <- cbind(1, x)
  C1 <- .fun(x, ...)
  M1 <- qr.Q(qr(cbind(C0, C1)))[, -(1:NCOL(C0))]
  S <- svd(t(M1) %*% C1 %*% M1)
  M2.inv.sqrt <- S$v %*% diag(1 / sqrt(S$d)) %*% t(S$u)
  K1 <- C1 %*% M1 %*% M2.inv.sqrt
  structure(
    list(C0 = C0, K1 = K1, M1 = M1, M2.inv.sqrt = M2.inv.sqrt),
    class = "ss.decomp"
    )
}





## basis.formals <- function(f) {
##   if (!is.function(f)) .Unrecognized("f", class(f))
##   if (!inherits(f, "dist.function")) {
##     if (length(b <- as.list(body(f))))
##       b <- if (b[[1L]] == as.name("{")) as.expression(b[-1L])
##         else as.expression(body(f))
##     else
##       b <- expression(NULL)
##     formals (f) <- alist(x = , y = , ... =)
##     body (f) <- as.call(c(as.name("{"), expression(if (missing(y)) y <- x), b))
##     class (f) <- c("dist.function", class(f))
##   }
##   f
## }










