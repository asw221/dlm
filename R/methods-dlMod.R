

#' @describeIn vcoef S4 Method for \code{"\link[=dlMod]{dlMod}"} Objects
setMethod("vcoef0", signature = "dlMod",
          function(object, scaled = TRUE, ...) {
            z <- c(lme4::getME(object, "beta"),
                   as.matrix(lme4::getME(object, "b"))
                   )
            if (scaled) c(as.matrix(scaleMat(object) %*% z)) else z
          })

#' @describeIn vcoef S4 Method for \code{"\link[=dlMod]{dlMod}"} Objects
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


#' @describeIn scaleMat S4 Method for \code{"\link[=dlMod]{dlMod}"} Objects
setMethod("scaleMat", signature = "dlMod",
          function(object, ...) {
            pq <- unlist(lme4::getME(object, c("p", "q")))
            S <- Matrix::Diagonal(sum(pq))
            ndx <- lagIndex(object)
            for (i in 1:length(ndx))
              S[ndx[[i]], ndx[[i]]] <- omega(object@bases[[object@index[i]]])
            S
          })




setMethod("cholfVar", signature = "dlMod",
  function(object, ...) {
    p <- lme4::getME(object, "p")
    q <- lme4::getME(object, "q")
    RZi <- t(solve(lme4::getME(object, "L"), Matrix::Diagonal(q), "L"))
    RXi <- solve(lme4::getME(object, "RX"))
    RZXi <- -RZi %*% lme4::getME(object, "RZX") %*% RXi
    Lambda <- lme4::getME(object, "Lambda")
    L <- Matrix::bdiag(RXi, Lambda %*% RZi)
    L[(p + 1):nrow(L), 1:p] <- Lambda %*% RZXi
    sigma(object) * L
  })



setMethod("Sigma", signature = "dlMod",
          function(object, scaled = TRUE, ...) {
            if (scaled)
              tcrossprod(scaleMat(object) %*% cholfVar(object))
            else
              tcrossprod(cholfVar(object))
          })


setMethod("changePoint", "dlMod",
          function(object, ...) {
            .Ignored(...)
            ci <- confint(object, ...)
            non0 <- !(ci[, 1] <= 0 & ci[, ncol(ci)] >= 0)
            lapply(lagIndex(object),
                   function(i) {
                     x <- non0[i]
                     which(x & !c(tail(x, -1), TRUE) & c(FALSE, head(x, -1)))
                   })
          })




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



