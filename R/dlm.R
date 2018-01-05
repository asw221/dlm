
## dlm
## -------------------------------------------------------------------
#' @title Distributed lag models
#'
#' @description
#' Fit distributed lag models using \pkg{lme4} to penalize smooth terms.
#' Other random effects terms and generalized mixed models supported.
#'
#' @usage
#' dlm(formula, data, subset, na.action, weights, offset,
#'     method = c("REML", "MLE"), family = gaussian(),
#'     control = list(), ...)
#'
#' @param formula
#'   an object of class \code{stats::\link[stats]{formula}}:
#'   a symbolic description of the model to be fitted. See Details
#' @param data
#'   an optional data frame, list, or environment containing the
#'   variables of the model to be fitted
#' @param subset
#'   optional vector specifying a subset of observations to be used in the
#'   fitting process
#' @param na.action
#'   optional function that indicates what should happen when the data contains
#'   \code{NA}'s. The default is set by the \code{na.action} setting of
#'   \code{base::\link[base]{options}}, usually
#'   \code{stats::\link[stats]{na.omit}}
#' @param weights
#'   optional vector of weights to be used in the fitting process. Should be
#'   \code{NULL} or a numeric vector
#' @param offset
#'   a known offset term to include in the model, e.g. for \code{poisson()}
#'   family models
#' @param method
#'   algorithm used to fit the DLM. Partial matching and capitalization allowed.
#'   The default is \code{"REML"} for linear/\code{gaussian(link = "identity")}
#'   family models, and \code{"MLE"} otherwise
#' @param family
#'   a description of the error distribution and link function to be used
#'   in the model. The default is \code{gaussian(link = "identity")}.
#'   See \code{stats::\link[stats]{family}} for possible family
#'   functions and details
#' @param control
#'   either a \code{list} object with arguments to be passed to the
#'   \code{lme4::\link[lme4]{lmerControl}} sequence, or the output of
#'   \code{[g]lmerConrol} directly
#' @param ...
#'   Additional parameters passed to \code{lme4::\link[lme4]{lFormula}}
#'
#' @details
#' Models are specified using typical \pkg{lme4} \code{formula} syntax with
#' at least one set of
#' lag terms returned by a given smoothing function (e.g. see \code{\link{cr}}).
#' The smoothing function can be any that returns a \code{\link{SmoothLag}} basis
#' object. See Examples
#' for a basic call to \code{dlm} using the formula interface, and a cubic
#' radial lag basis specified via \code{cr}, and the \pkg{dlmBE}
#' \code{\link[=dlmBE-package]{package documentation}} for a discussion
#' of the types of models \code{dlm} is designed to handle.
#'
#'
#' @return
#' An S4 object that inherits from \code{\link{dlMod}} (and
#' \code{lme4::\link[lme4]{merMod}}, by extension) containing the
#' results of the fitted model. Many standard model summary methods are
#' available for these object types
#'
#'
#' @references Baek J, Sanchez BN, Berrocal VJ, & Sanchez-Vaznaugh EV
#' (2016) Epidemiology 27(1):116-24.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26414942}{PubMed})
#'
#' @references Baek J, Hirsch JA, Moore K, Tabb LP, et al. (2017)
#' Epidemiology 28(3):403-11.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/28145983}{PubMed})
#'
#' @references Bates D, Maechler M, Bolker BM, & Walker SC (2015)
#' Fitting linear mixed-effects models using lme4. J
#' Stat Softw 67(1).
#' (\href{https://www.jstatsoft.org/article/view/v067i01/0}{jstatsoft.org})
#'
#'
#' @examples
#' data (simdata)
#'
#' ## Setup distance count matrix and corresponding lag distances
#' X <- as.matrix(simdata[, -(1:3)])
#' lag <- seq(0.1, 10, length.out = ncol(X))
#'
#' fit <- dlm(Y ~ Age + Gender + cr(lag, X), data = simdata)
#' summary (fit)
#'
#' @seealso \code{lme4::\link[lme4]{lmer}}, \code{\link{cr}},
#'   \code{\link{dlMod}}
#'
dlm <- function(formula, data, subset, na.action, weights, offset,
  method = c("REML", "MLE", "reml", "mle"),
  family = gaussian(),
  control = list(),
  ...
) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"),
             names(mf), 0
             )
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- lme4::subbars(formula)
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  parsed <- interpret.dlm(formula, mf)
  rm (mf, m)

  method <- tolower(match.arg(method))
  if (method == "reml" || method == "mle") {
    reml <- method == "reml"
    fit <- lme4.dlm(parsed, family, control, REML = reml, ...)
  }
  else
    stop ("Unrecognized method, ", method)
  makeDlMod(fit, parsed, match.call())
}
## dlm




## lme4.dlm
## -------------------------------------------------------------------
#' @title Fit distributed lag models using lme4
#'
#' @description
#' Fits an interpreted distributed lag model using \pkg{lme4}
#' \code{\link[lme4]{modular}} functions
#'
#' @param parsed
#'   an interpreted \code{dlm} formula object returned by
#'   \code{\link{interpret.dlm}}
#' @param family
#'   a description of the error distribution and link function to be used
#'   in the model. The default is \code{gaussian(link = "identity")}.
#'   See \code{stats::\link[stats]{family}} for possible family
#'   functions and details
#' @param control
#'   either a \code{list} object with arguments to be passed to the
#'   \code{lme4::\link[lme4]{lmerControl}} sequence, or the output of
#'   \code{[g]lmerConrol} directly
#' @param REML
#'   if \code{TRUE} and a linear \code{\link{dlm}} model is specified,
#'   \code{lme4.dlm} will use REML to fit the model. MLE will be
#'   use otherwise
#' @param ...
#'   other parameters to be passed to the \pkg{lme4}
#'   \code{\link[lme4]{modular}} family functions
#'
#' @details
#' Together with \code{\link{interpret.dlm}}, this function does the
#' main grunt work for \code{\link{dlm}}. Given an interpreted model,
#' \code{lme4.dlm} organizes the parsed data into the \pkg{lme4}
#' \code{\link[lme4]{modular}} functions to fit the model and returns
#' the fit as an \pkg{lme4} object.
#'
#' @return an object that inherits from \code{lme4::\link[lme4]{merMod}}
#'   containing a fitted model
#'
lme4.dlm <- function(parsed, family = gaussian(),
                     control = list(), REML = FALSE,
                     ...
) {
  if (!inherits(parsed, "parsed.dlm")) .Unrecognized("parsed", class(parsed))
  if (is.function(family))
    family <- family()
  linear <- (family$family == "gaussian" && family$link == "identity")

  if (!inherits(control, "merControl")) {
    if (is.list(control)) {
      control <- if (linear) do.call(lme4::lmerControl, control)
                 else do.call(lme4::glmerControl, control)
    }
    else
      .Unrecognized("control", class(control))
  }

  ## reassign environment for formula so [g]lFormula knows where
  ## to look for variable parsed in calls to eval()
  ##   -- may be a better way to do this!
  formula <- parsed$lme4.formula
  environment (formula) <- sys.frame(sys.nframe())

  ## Figure out the correct lme4 modular functions to fit the model
  if (linear) {
    ## weights and offset are NULL values if missing and are handled
    ## correctly in any case
    pf <- lme4::lFormula(formula, data = parsed$model,
                   control = control, REML = REML,
                   weights = parsed$model[["(weights)"]],
                   offset = parsed$model[["(offset)"]],
                   ...
                   )
    .Deviance <- lme4::mkLmerDevfun
    .Optimize <- lme4::optimizeLmer
  }
  else {  # glFormula needs a 'family' argument; all else the same
    ## REML only for linear case
    pf <- lme4::glFormula(formula, data = parsed$model,
                   control = control, family = family,
                   weights = parsed$model[["(weights)"]],
                   offset = parsed$model[["(offset)"]],
                   ...
                   )
    .Deviance <- lme4::mkGlmerDevfun
    .Optimize <- lme4::optimizeGlmer
  }

  ## Replace iid dummy variables corresponding to penalized lag terms
  ## in RE design matrix with the actual lag term data
  lnms <- attr(parsed$lag.group, "dictionary")
  m <- match(names(lnms), names(pf$reTrms$cnms))
  pf$reTrms <- within(pf$reTrms, {
    for (i in seq_along(lnms)) {
      j <- m[i]
      Zt[(Gp[j]+1):Gp[j+1], ] <- parsed$Bt[parsed$lag.group[lnms[i]], ]
      names (cnms)[j] <- lnms[i]
      cnms[[i]] <- "(mean)"
    }
  })
  ## pf$reTrms <- within(pf$reTrms, {
  ##   for (i in seq_along(lnms)) {
  ##     Zt[(Gp[i]+1):Gp[i+1], ] <- parsed$Bt[parsed$lag.group[lnms[m[i]]], ]
  ##     names (cnms)[i] <- lnms[m[i]]
  ##     cnms[[i]] <- "(mean)"
  ##   }
  ## })
  devfun <- do.call(.Deviance, pf)
  optim <- .Optimize(devfun)
  fit <- lme4::mkMerMod(rho = environment(devfun),
                  opt = optim, reTrms = pf$reTrms, fr = pf$fr
                  )
  names (fit@flist) <- names(fit@cnms)
  attr(fit@frame, "formula") <- parsed$lme4.formula
  fit
}
## lme4.dlm




## makeDlMod.merMod
## -------------------------------------------------------------------
#' @title Convert fitted models to 'dlMod' objects
#'
#' @description
#' Convert an appropriately fit model object into an object of class
#' \code{\link{dlMod}}
#'
#' @param object
#'   a fitted model object
#' @param parsed
#'   an interpreted \code{dlm} formula object returned by
#'   \code{\link{interpret.dlm}}
#' @param call
#'   an optional matched function call
#'
#' @return an S4 object of class \code{\link{dlMod}}
#'
#' @name makeDlMod
makeDlMod.merMod <- function(object, parsed, call, ...) {
  if (!(lme4::isLMM(object) || lme4::isGLMM(object)))
    stop ("Not yet implemented for ", class(object), " objects")
  if (!inherits(parsed, "parsed.dlm")) .Unrecognized("parsed", class(parsed))
  if (missing(call)) call <- object@call
  else if (!is.call(call)) .Unrecognized("call", class(call))
  index <- parsed$bi
  names (index) <- attr(parsed$lag.group, "dictionary")
  obj <- dlMod(object, bases = parsed$bases, index = index, ...)
  obj@call <- call
  obj@resp <- object@resp
  obj
}
## makeDlMod.merMod




