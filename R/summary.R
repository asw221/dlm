


## from lme4:
summary.dlMod <- function (object,
    correlation = (p <= getOption("lme4.summary.cor.max")),
    use.hessian = NULL, ...)
{
  ## Most of this is lifted straight from summary.merMod()
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
  structure(list(methTitle = methTitle(dd), objClass = class(object),
                 devcomp = devC, isLmer = is(resp, "lmerResp"), useScale = useSc,
                 logLik = llAIC[["logLik"]], family = famL$family, link = famL$link,
                 ngrps = ngrps(object), coefficients = coefs, sigma = sig,
                 vcov = vcov(object, correlation = correlation, sigm = sig),
                 varcor = varcor, AICtab = llAIC[["AICtab"]], call = object@call,
                 residuals = residuals(object, "pearson", scaled = TRUE),
                 fitMsgs = msgs(object), optinfo = object@optinfo),
            class = c("summary.dlMod", "summary.merMod"))
}




print.summary.merMod <- function (
    x, digits = max(3, getOption("digits") - 3), correlation = NULL,
    symbolic.cor = FALSE, signif.stars = getOption("show.signif.stars"),
    ranef.comp = c("Variance", "Std.Dev."), show.resids = TRUE,
    ...)
{
    .prt.methTit(x$methTitle, x$objClass)
    .prt.family(x)
    .prt.call(x$call)
    cat("\n")
    .prt.aictab(x$AICtab)
    cat("\n")
    if (show.resids)
        .prt.resids(x$residuals, digits = digits)
    .prt.VC(x$varcor, digits = digits, useScale = x$useScale,
        comp = ranef.comp, ...)
    .prt.grps(x$ngrps, nobs = x$devcomp$dims[["n"]])
    p <- nrow(x$coefficients)
    if (p > 0) {
        cat("\nFixed effects:\n")
        printCoefmat(x$coefficients, zap.ind = 3, digits = digits,
            signif.stars = signif.stars)
        hasCor <- !is.null(VC <- x$vcov) && !is.null(VC@factors$correlation)
        if (is.null(correlation)) {
            cor.max <- getOption("lme4.summary.cor.max")
            correlation <- hasCor && p <= cor.max
            if (!correlation && p > cor.max) {
                nam <- deparse(substitute(x))
                if (length(nam) > 1 || nchar(nam) >= 32)
                  nam <- "...."
                message(sprintf(paste("\nCorrelation matrix not shown by default, as p = %d > %d.",
                  "Use print(%s, correlation=TRUE)  or", "\t vcov(%s)\t if you need it\n",
                  sep = "\n"), p, cor.max, nam, nam))
            }
        }
        else if (!is.logical(correlation))
            stop("'correlation' must be NULL or logical")
        if (correlation) {
            if (is.null(VC))
                VC <- vcov(x, correlation = TRUE)
            corF <- VC@factors$correlation
            if (is.null(corF)) {
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
                  dimnames(corf) <- list(rns, abbreviate(rn,
                    minlength = 1, strict = TRUE))
                  print(symnum(corf))
                }
                else {
                  corf <- matrix(format(round(corF@x, 3), nsmall = 3),
                    ncol = p, dimnames = list(rns, abbreviate(rn,
                      minlength = 6)))
                  corf[!lower.tri(corf)] <- ""
                  print(corf[-1, -p, drop = FALSE], quote = FALSE)
                }
            }
        }
    }
    if (length(x$fitMsgs) && any(nchar(x$fitMsgs) > 0)) {
        cat("fit warnings:\n")
        writeLines(x$fitMsgs)
    }
    .prt.warn(x$optinfo, summary = FALSE)
    invisible(x)
}

