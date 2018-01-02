

## These functions are lifted almost entirely straight from lme4
## Minor modifications only for dlMod objects:
##   - remove printing of "fixed" part of lag terms in
##     summary.dlMod coefficient matrix
##

summary.dlMod <- function (object,
    correlation = (p <= getOption("lme4.summary.cor.max")),
    use.hessian = NULL, ...)
{
  if (length(list(...)) > 0) {
    warning ("additional arguments ignored")
  }
  msgs <- function(x) { # use on object@pp$X
    aX <- attributes(x)
    .msgs <- unlist(aX[grepl("^msg", names(aX))])
    if (!is.null(.msgs)) .msgs else character()
  }
  hess.avail <- (!is.null(h <- object@optinfo$derivs$Hessian) &&
                 nrow(h) > length(lme4::getME(object, "theta")))
  if (is.null(use.hessian))
    use.hessian <- hess.avail
  if (use.hessian && !hess.avail)
    stop("'use.hessian=TRUE' specified, but Hessian is unavailable")
  resp <- object@resp
  devC <- object@devcomp
  dd <- devC$dims
  useSc <- as.logical(dd[["useSc"]])
  sig <- sigma(object)
  famL <- if (is(resp, "glmResp")) resp$family[c("family", "link")]
    else list(family = NULL, link = NULL)
  p <- length(coefs <- fixef(object))
  vc <- vcov(object, use.hessian = use.hessian)
  stdError <- sqrt(diag(vc))
  coefs <- cbind(Estimate = coefs, `Std. Error` = stdError)
  if (p > 0) {
    coefs <- cbind(coefs, (cf3 <- coefs[, 1]/coefs[, 2]),
                   deparse.level = 0)
    colnames(coefs)[3] <- paste(if (useSc)
                                  "t"
                                else "z", "value")
    if (isGLMM(object))
      coefs <- cbind(coefs, `Pr(>|z|)` = 2 * pnorm(abs(cf3),
                                                   lower.tail = FALSE))
  }
  llAIC <- lme4::llikAIC(object)
  varcor <- lme4::VarCorr(object)
  structure(list(methTitle = methTitle(dd),
                 objClass = class(object),
                 devcomp = devC,
                 isLmer = is(resp, "lmerResp"),
                 useScale = useSc,
                 logLik = llAIC[["logLik"]],
                 family = famL$family,
                 link = famL$link,
                 ngrps = ngrps(object),
                 coefficients = coefs,
                 sigma = sig,
                 vcov = vcov(object, correlation = correlation, sigm = sig),
                 varcor = varcor,
                 AICtab = llAIC[["AICtab"]],
                 call = object@call,
                 residuals = residuals(object, "pearson", scaled = TRUE),
                 fitMsgs = msgs(object),
                 optinfo = object@optinfo,
                 lag.names = names(object@index)
                 ),
            class = c("summary.dlMod", "summary.merMod"))
}
## summry.dlMod



print.summary.dlMod <- function(x,
  digits = max(3, getOption("digits") - 3),
  correlation = NULL,
  symbolic.cor = FALSE,
  signif.stars = getOption("show.signif.stars"),
  ranef.comp = c("Variance", "Std.Dev."),
  show.resids = TRUE,
  ...
) {
  .prt.methTit(x$methTitle, x$objClass)
  .prt.family(x)
  .prt.call(x$call); cat("\n")
  .prt.aictab(x$AICtab); cat("\n")
  if (show.resids)
    ## need residuals.merMod() rather than residuals():
    ##  summary.merMod has no residuals method
    .prt.resids(x$residuals, digits = digits)
  .prt.VC(x$varcor, digits = digits, useScale = x$useScale,
          comp = ranef.comp, ...)
  .prt.grps(x$ngrps, nobs = x$devcomp$dims[["n"]])

  p <- nrow(x$coefficients)
  ## find index of non-lag covariates (ic)
  ic <- parse.names(x$lag.names, rownames(x$coefficients), .warn = FALSE)[""]
  if (p > 0 && length(ic) > 0) {
    cat("\nFixed effects:\n")
    printCoefmat(x$coefficients[ic, , drop = FALSE], zap.ind = 3, #, tst.ind = 4
                 digits = digits, signif.stars = signif.stars)
    ## do not show correlation when   summary(*, correlation=FALSE)  was used:
    hasCor <- !is.null(VC <- x$vcov) && !is.null(VC@factors$correlation)
    if(is.null(correlation)) { # default
      cor.max <- getOption("lme4.summary.cor.max")
      correlation <- hasCor && p <= cor.max
      if(!correlation && p > cor.max) {
        nam <- deparse(substitute(x))
        if(length(nam) > 1 || nchar(nam) >= 32) nam <- "...."
        message(sprintf(paste(
          "\nCorrelation matrix not shown by default, as p = %d > %d.",
          "Use print(%s, correlation=TRUE)  or",
          "	 vcov(%s)	 if you need it\n", sep = "\n"),
          p, cor.max, nam, nam))
      }
    }
    else if(!is.logical(correlation)) stop("'correlation' must be NULL or logical")
    if(correlation) {
      if(is.null(VC)) VC <- vcov(x, correlation = TRUE)
      corF <- VC@factors$correlation
      if (is.null(corF)) { # can this still happen?
        message("\nCorrelation of fixed effects could have been required in summary()")
        corF <- cov2cor(VC)
      }
      p <- ncol(corF)
      if (p > 1) {
        rn <- rownames(x$coefficients)
        rns <- abbreviate(rn, minlength = 11)
        cat("\nCorrelation of Fixed Effects:\n")
        if (is.logical(symbolic.cor) && symbolic.cor) {
          corf <- as(corF, "matrix")
          dimnames(corf) <- list(rns,
                                 abbreviate(rn, minlength = 1, strict = TRUE))
          print(symnum(corf))
        } else {
          corf <- matrix(format(round(corF@x, 3), nsmall = 3),
                         ncol = p,
                         dimnames = list(rns, abbreviate(rn, minlength = 6)))
          corf <- corf[ic, ic, drop = FALSE]
          corf[!lower.tri(corf)] <- ""
          print(corf[-1, -NCOL(corf), drop = FALSE], quote = FALSE)
        } ## !symbolic.cor
      }  ## if (p > 1)
    } ## if (correlation)
  } ## if (p>0)

  if(length(x$fitMsgs) && any(nchar(x$fitMsgs) > 0)) {
    cat("fit warnings:\n"); writeLines(x$fitMsgs)
  }
  .prt.warn(x$optinfo,summary=FALSE)
  invisible(x)
}
## print.summary.dlMod



## Formatting utilities (also from lme4)
## -------------------------------------------------------------------


cat.f <- function(...) cat(..., fill = TRUE)

.prt.methTit <- function(mtit, class) {
  if(nchar(mtit) + 5 + nchar(class) > (w <- getOption("width"))) {
    ## wrap around
    mtit <- strwrap(mtit, width = w - 2, exdent = 2)
    cat(mtit, " [",class,"]", sep = "", fill = TRUE)
  } else ## previous: simple one-liner
    cat(sprintf("%s ['%s']\n", mtit, class))
}

.prt.family <- function(famL) {
  if (!is.null(f <- famL$family)) {
    cat.f(" Family:", f,
          if(!is.null(ll <- famL$link)) paste(" (", ll, ")"))
  }
}

.prt.resids <- function(resids, digits, title = "Scaled residuals:", ...) {
  cat(title,"\n")
  ## FIXME: need testing code
  rq <- setNames(zapsmall(quantile(resids, na.rm=TRUE), digits + 1L),
                 c("Min", "1Q", "Median", "3Q", "Max"))
  print(rq, digits = digits, ...)
  cat("\n")
}

.prt.call <- function(call, long = TRUE) {
  if (!is.null(cc <- call$formula))
    cat.f("Formula:", deparse(cc))
  if (!is.null(cc <- call$data))
    cat.f("   Data:", deparse(cc))
  if (!is.null(cc <- call$weights))
    cat.f("Weights:", deparse(cc))
  if (!is.null(cc <- call$offset))
    cat.f(" Offset:", deparse(cc))
  if (long && length(cc <- call$control) &&
      !identical((dc <- deparse(cc)), "lmerControl()"))
    ## && !identical(eval(cc), lmerControl()))
    cat.f("Control:", dc)
  if (!is.null(cc <- call$subset))
    cat.f(" Subset:", deparse(cc))
}


.prt.aictab <- function(aictab, digits = 1) {
  t.4 <- round(aictab, digits)
  if (length(aictab) == 1 && names(aictab) == "REML")
    cat.f("REML criterion at convergence:", t.4)
  else {
    ## slight hack to get residual df formatted as an integer
    t.4F <- format(t.4)
    t.4F["df.resid"] <- format(t.4["df.resid"])
    print(t.4F, quote = FALSE)
  }
}

.prt.VC <- function(varcor, digits, comp, formatter = format, ...) {
  cat("Random effects:\n")
  fVC <- if(missing(comp))
           formatVC(varcor, digits = digits, formatter = formatter)
         else
           formatVC(varcor, digits = digits, formatter = formatter, comp = comp)
  print(fVC, quote = FALSE, digits = digits, ...)
}



.prt.grps <- function(ngrps, nobs) {
  cat(sprintf("Number of obs: %d, groups: ", nobs),
      paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "),
      fill = TRUE)
}



.prt.warn <- function(optinfo, summary=FALSE, ...) {
  if(length(optinfo) == 0) return() # currently, e.g., from refitML()
  ## check all warning slots: print numbers of warnings (if any)
  cc <- optinfo$conv$opt
  msgs <- unlist(optinfo$conv$lme4$messages)
  ## can't put nmsgs/nwarnings compactly into || expression
  ##   because of short-circuiting
  nmsgs <- length(msgs)
  warnings <- optinfo$warnings
  nwarnings <- length(warnings)
  if (cc > 0 || nmsgs > 0 || nwarnings > 0) {
    if (summary) {
      cat(sprintf("convergence code %d; %d optimizer warnings; %d lme4 warnings",
                  cc,nmsgs,nwarnings),"\n")
    } else {
      cat(sprintf("convergence code: %d", cc),
          msgs,
          unlist(warnings),
          sep="\n")
      cat("\n")
    }
  }
}



