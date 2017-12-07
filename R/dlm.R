
#' @title Distributed lag models
#'
#' @description
#' Fit distributed lag models
#'
#' @usage
#' dlm(formula, data, subset, na.action, weights, offset,
#'     method = c("REML", "MLE"), family = gaussian(),
#'     control = list(), ...)
#'
#' @param formula
#'   an object of class \code{"\link[stats]{formula}"}:
#'   a symbolic description of the model to be fitted. See Details.
#' @param data
#'   an optional data frame, list, or environment containing the
#'   variables of the fitted model.
#' @param subset
#'   optional vector specifying a subset of observations to be used in the
#'   fitting process.
#' @param na.action
#'   optional function that indicates what should happen when the data contains
#'   \code{NA}'s.
#' @param method
#'   method used to fit the DLM. Partial matching and capitalization allowed.
#' @param control
#'   a list of simulation control parameters. See Details.
#' @param ...
#'   Additional parameters passed to \code{\link[lme4]{lFormula}}.
#'
#' @details
#' Models are specified using familiar R \code{formula} syntax with one set of
#' lag terms returned by a given smoothing function (e.g. see \code{\link{cr}}).
#' The smoothing function can be any that returns a \code{\link{SmoothLag}} basis
#' object. Multiple lag terms or interactions with lag terms are not allowed, nor
#' are \code{NA} and missing-value lag terms supported. See Examples
#' for a basic call to \code{dlm} using the formula interface, and a cubic
#' radial lag basis specified via \code{cr}.
#'
#' The \code{control} list specifies additional optional \code{"bayes"} method
#' arguments, and may include: \code{n.sims}, the total number of simulations
#' to run (default = 5000); \code{n.save}, the total number of simulations to
#' save (default = 1000); \code{burnin}, number of simulations to discard from
#' the start of the chain (default = 2000);
#' \code{alpha.tau.prior}, prior (shape) parameter
#'   \eqn{\alpha_{\tau^2}}{\alpha[\tau^2]} (default = 0.1);
#' \code{beta.tau.prior}, prior (rate) parameter
#'   \eqn{\beta_{\tau^2}}{\alpha[\tau^2]} (default = 1e-6);
#' \code{alpha.sigma.prior}, prior (shape) parameter
#'   \eqn{\alpha_{\sigma^2}}{\alpha[\sigma^2]} (default = 0.1);
#' \code{beta.sigma.prior}, prior (rate) parameter
#'   \eqn{\beta_{\sigma^2}}{\beta[\sigma^2]} (default = 1e-6).
#'
#' The prior distribution heirarchy we assume in the Bayesian
#' formulation of the DLM is as follows:
#'
#' \deqn{y \sim N(D \theta, \sigma^2 I_n)}{y ~ N(D * \theta, \sigma^2 * In)}
#' \deqn{\theta \sim N(\mu_{\theta}, \Sigma_{\theta})}{\theta ~ N(\mu, \Sigma)}
#' \deqn{\sigma^2 \sim \mathrm{Inv-Gamma}(\alpha_{\sigma^2}, \frac{1}{\beta_{\sigma^2}})}{\sigma^2 ~ Inv-Gamma(\alpha[\sigma^2], 1 / \beta[\sigma^2])}
#'
#' with,
#' \deqn{\theta_l \sim N(\mu_l, \tau^2)}{\theta(l) ~ N(\mu(l), \tau^2)}
#' \deqn{\tau^2 \sim \mathrm{Inv-Gamma}(\alpha_{\tau^2}, \frac{1}{\beta_{\tau^2}})}{\tau^2 ~ Inv-Gamma(\alpha[\tau^2], 1 / \beta[\tau^2])}
#'
#' where \eqn{l \in L}{l :- L} indexes the set of lag coefficients.
#'
#' @return
#' An S4 object that inherits from \code{"\link[=Dlm-class]{Dlm}"} containing
#' the results of the fitted model. If construction of this object fails,
#' \code{dlm} will issue a warning, and as a last resort attempt to return a
#' list with components of the fitted model.
#'
#'
#' @references Baek J, Sanchez BN, Berrocal BJ, & Sanchez-Vaznaugh EV
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
#' @seealso \code{\link[lme4]{lmer}}, \code{\link{cr}},
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
      Zt[(Gp[i]+1):Gp[i+1], ] <- parsed$Bt[parsed$lag.group[lnms[m[i]]], ]
      names (cnms)[i] <- lnms[m[i]]
      cnms[[i]] <- "(mean)"
    }
  })
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




