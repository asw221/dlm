
## Sigma
## -------------------------------------------------------------------
#' @rdname estimands
setMethod("Sigma", signature = "dlMod",
          function(object, scaled = TRUE, ...) {
            if (scaled)
              tcrossprod(scaleMat(object) %*% cholfVar(object))
            else
              tcrossprod(cholfVar(object))
          })
## Sigma

## vcoef
## -------------------------------------------------------------------
#' @rdname estimands
setMethod("vcoef0", signature = "dlMod",
          function(object, scaled = TRUE, ...) {
            z <- c(lme4::getME(object, "beta"),
                   as.matrix(lme4::getME(object, "b"))
                   )
            if (scaled) c(as.matrix(scaleMat(object) %*% z)) else z
          })
## vcoef0

#' @rdname estimands
setMethod("vcoef", signature = "dlMod",
          function(object, scaled = TRUE, ...) {
            z <- vcoef0(object, scaled = scaled)
            ## build names(z):
            p <- lme4::getME(object, "p")
            ng <- diff(object@Gp)
            nms <- character(length(z))
            nms[1:p] <- colnames(lme4::getME(object, "X"))
            for (i in seq_along(object@cnms)) {
              nm <- names(object@cnms)[i]
              if (nm %in% names(object@index)) {  # spline term
                ## paste in decomposed lag term scales
                ## the first few are treated as fixed effects
                x <- tail(object@bases[[object@index[nm]]]@x, ng[i])
                new.nms <- paste(nm, x, sep = "")
              }
              else {  # non-spline random effect
                ## paste in factor levels
                new.nms <- paste(nm, levels(object@flist[[i]]), sep = "")
                new.nms <- sapply(new.nms, paste, object@cnms[[i]], sep = ".")
              }
              nms[(object@Gp[i] + 1):object@Gp[i + 1] + p] <- new.nms
            }
            structure(z, names = nms)
          })
## vcoef


## scaleMat
## -------------------------------------------------------------------
#' @title Extract distributed lag scale matrix
#'
#' @description
#' Return lag coefficient scale matrix, \eqn{S}. \eqn{S} should be invertable.
#' Typically users should only interact with \code{scaleMat} indirectly
#' through the \code{\link{estimands}} functions (with the argument
#' \code{scaled = TRUE}).
#'
#' @param object a fitted model object
#'
#' @details
#' If \eqn{\theta} is the vector of unscaled regression coefficients,
#' then corresponding distributed lag terms in the vector
#' \eqn{S \theta}{S * \theta} can be interpreted as average regression
#' coefficients up to their associated radius. Non-DL coefficients are
#' unchanged by this transformation.
#'
#' @return A square numeric matrix represented as a
#' \code{Matrix::\link[Matrix]{sparseMatrix}} object
#'
#' @name scaleMat
setMethod("scaleMat", signature = "dlMod",
          function(object, ...) {
            pq <- unlist(lme4::getME(object, c("p", "q")))
            S <- Matrix::Diagonal(sum(pq))
            ndx <- lagIndex(object)
            for (i in 1:length(ndx))
              S[ndx[[i]], ndx[[i]]] <- omega(object@bases[[object@index[i]]])
            S
          })



## changePoint
## -------------------------------------------------------------------
#' @title Lag coefficient change points
#'
#' @description
#' For each set of distributed lag terms in the model, finds and returns
#' the (named) integer indices of radii where the corresponding
#' coefficient is significantly different from zero, but the immediately
#' larger radius is not significantly different from zero.
#'
#' @param object
#'   a fitted model object with viable \code{\link{lagIndex}} and
#'   \code{\link[=estimands]{confint}} methods
#' @param ...
#'   additional arguments passed to \code{\link[=estimands]{confint}}
#'
#'
#' @return A list of integer vectors. One list element for each set of
#'   DL terms in the model
#'
#' @examples
#' data (simdata)
#'
#' ## Setup distance count matrix and corresponding lag distances
#' X <- as.matrix(simdata[, -(1:3)])
#' lag <- seq(0.1, 10, length.out = ncol(X))
#'
#' fit <- dlm(Y ~ Age + Gender + cr(lag, X), data = simdata)
#' changePoint(fit)
#'
#' @name changePoint
setMethod("changePoint", "dlMod",
          function(object, ...) {
            ci <- confint(object, coef = FALSE, ...)
            non0 <- !(ci[, 1] <= 0 & ci[, ncol(ci)] >= 0)
            lapply(lagIndex(object),
                   function(i) {
                     x <- non0[i]
                     which(x & !c(tail(x, -1), TRUE) & c(FALSE, head(x, -1)))
                   })
          })
## changePoint


## lagIndex
## -------------------------------------------------------------------
#' @title Extract list of indices of lag terms
#'
#' @description
#' Intended for use in conjunction with the \code{\link{estimands}} functions.
#' For each set of distributed lag terms in the model, find and return the
#' integer indices that would extract the corresponding coefficients from
#' the vector returned by \code{vcoef}.
#'
#' @param object
#'   a fitted model object
#' @param .fixed
#'   a logical flag (default \code{.fixed = TRUE}) to indicate whether the
#'   fixed/unpenalized DL term coefficient indices should be included or not.
#'
#' @examples
#' data (simdata)
#'
#' ## Setup distance count matrix and corresponding lag distances
#' X <- as.matrix(simdata[, -(1:3)])
#' lag <- seq(0.1, 10, length.out = ncol(X))
#'
#' fit <- dlm(Y ~ Age + Gender + cr(lag, X), data = simdata)
#' vcoef(fit)[lagIndex(fit)[[1]]]
#'
#' @seealso \code{\link{vcoef}}
#' @name lagIndex
setMethod("lagIndex", "dlMod",
          function(object, .fixed = TRUE, ...) {
            .ind <- function(obj, name) {
              i <- which(names(obj@cnms) == name)
              (obj@Gp[i] + 1):obj@Gp[i + 1] + lme4::getME(obj, "p")
            }
            .Ignored(...)
            lnms <- names (object@index)
            ndx <- lapply(lnms, function(nm) .ind(object, nm))
            names (ndx) <- lnms
            if (.fixed) {
              fe.nms <- colnames(lme4::getME(object, "X"))
              grp <- parse.names(lnms, fe.nms, .warn = FALSE)
              for (nm in lnms)
                ndx[[nm]] <- c(grp[nm], ndx[[nm]])
            }
            ndx
          })
## lagIndex



## coef.dlMod
## -------------------------------------------------------------------
#' @rdname estimands
coef.dlMod <- function(object, scaled = TRUE, ...) {
  .Ignored(...)
  z <- vcoef(object, scaled = scaled)
  if (length(w <- which(!(names(object@cnms) %in% names(object@index))))) {
    p <- lme4::getME(object, "p")
    li <- lagIndex(object, .fixed = FALSE)
    fef <- t(z[c(1:p, unlist(li))])
    val <- list(length(w))
    names (val) <- names(object@cnms[w])
    for (i in seq_along(w)) {
      j <- w[i]
      refi <- matrix(z[(object@Gp[j] + 1):object@Gp[j + 1] + p],
                     ncol = length(object@cnms[[j]]), byrow = TRUE)
      colnames (refi) <- object@cnms[[j]]
      val[[i]] <- data.frame(rep(1, nrow(refi)) %*% fef, check.names = FALSE)
      rownames (val[[i]]) <- as.character(levels(object@flist[[j]]))
      for (nm in colnames(refi))
        val[[i]][[nm]] <- if (is.null(val[[i]][[nm]])) refi[, nm]
          else val[[i]][[nm]] + refi[, nm]
    }
    z <- val
  }
  structure(z, class = "coef.mer")
}
## coef.dlMod


## confint.dlMod
## -------------------------------------------------------------------
#' @param parm
#'   an integer or character index to subset parameters
#' @param level
#'   the desired confidence level
#' @param coef
#'   if \code{coef = TRUE} (the default), \code{confint.dlMod} will
#'   include an extra column (with \code{colname} \code{"coef"}) for the
#'   regression coefficients themselves
#'
#' @rdname estimands
confint.dlMod <- function(object, parm, level = 0.95, scaled = TRUE,
                          coef = TRUE, ...) {
  .Ignored(...)
  if (any(level >= 1 | level <= 0))  stop ("level out of bounds")
  b <- vcoef(object, scaled = scaled)
  a <- (1 - level) / 2
  a <- sort(c(a, 1 - a))
  q <- qnorm(a)
  se <- sqrt(diag(Sigma(object, scaled = scaled)))
  ci <- b + se %o% q
  colnames (ci) <- sprintf("%.1f%%", a * 100)
  rownames (ci) <- names(b)
  if (coef)
    ci <- cbind(coef = b, ci)
  if (!missing(parm)) ci[parm, ] else ci
}
## confint.dlMod



