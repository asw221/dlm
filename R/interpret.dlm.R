
## interpret.dlm
## -------------------------------------------------------------------
#' @title Interpret a DLM formula
#'
#' @description
#' Given an appropriate DL model formula and data-frame, prepares
#' the data to be fit by interpreting the formula and extracting
#' available smooth lag terms.
#'
#' @param formula
#'   a symbolic description of the model to be fitted.
#'   See \code{\link{dlm}} for details.
#' @param data
#'   a model-frame containing the data for each term in the model.
#'   Should already be appropriately subset, etc.
#' @param .names.func
#'   a function for creating names of dummy variables that act as
#'   placeholders for penalized spline terms in \pkg{lme4}'s setup.
#'   There should not be a need to alter this in normal use cases
#'
#' @details
#' Users should not typically have to interact with \code{interpret.dlm}
#' directly, but it may be useful for extensions.
#'
#' Uses \code{R}'s \code{stats::\link[stats]{model.matrix}} mechanisms to build
#' and parse the random effects (or penalized) components of spline-lag terms
#' in the model. The object returned is later passed to other \pkg{dlmBE}
#' functions in order to fit the specified model.
#'
#'
#' @return
#' an S3 object of class \code{"parsed.dlm"} with list elements:
#' \describe{
#'   \item{\code{formula}}{the formula passed to \code{interpret.dlm}}
#'   \item{\code{lme4.formula}}{a reconstructed formula that is then passed
#'     to the \pkg{lme4} \code{\link[lme4]{modular}} functions}
#'   \item{\code{model}}{a \code{data.frame} returned by call to
#'     \code{stats::\link[stats]{model.frame}}; a copy of the fixed
#'     effects data needed to fit the model}
#'   \item{\code{Bt}}{a matrix of the random or penalized lag basis vectors,
#'     where each vector is a row. Stored as an object that inherits from
#'     \code{Matrix::\link[Matrix]{dMatrix}}}
#'   \item{\code{bases}}{a list of all the unique bases represented in the
#'     \code{formula}. This may be \eqn{\leq}{<=} the number of separate
#'     spline-lag terms. All elements should inherit from
#'     \code{\link{SmoothLag}}}
#'   \item{\code{lag.group}}{an integer vector returned by
#'     \code{\link{parse.names}} where each unique integer corresponds to a
#'     separate spline-lag term. For lag term \code{i}, \code{lag.group == i}
#'     indexes the rows of \code{Bt} that correspond to the set of random or
#'     penalized basis vectors for that term}
#'   \item{\code{bi}}{for "basis index." Each set of lag terms indexed in
#'     \code{lag.group} has a matching basis decomposition in \code{bases}.
#'     \code{bi} keeps track of that matching}
#' }
#'
interpret.dlm <- function(formula, data,
  .names.func = function(n) paste("pseudoGroups", n, sep = "")
) {
  fake.formula <- Reduce(paste, deparse(formula))
  mt <- terms(data)
  mtf <- attr(mt, "factor")
  ## attr(mt.tmp, "intercept") <- 0  ## RE design does not have intercept
  if (is.empty.model(mt))
    stop ("Empty model")

  ## to use R's own model.matrix() mechanisms for building main
  ## effects and interactions with DL terms, duplicate the model-frame
  ## but replace the "main effects" components of smooth DL terms with
  ## their "random effects" basis sets
  mf.lag <- data
  lag <- logical(NROW(mtf))  ## which terms in data are lag terms
  bases <- list()
  for (j in 1:NROW(mtf))  if (lag[j] <- inherits(data[[j]], "SmoothLag")) {
    mf.lag[[j]] <- mf.lag[[j]]@random[attr(data, "na.action"), ]
    bases[[j]] <- data[[j]]@basis
    data[[j]] <- data[[j]]@.Data
  }
  if (!any(lag))
    stop ("Model has no apparent DL terms")

  bases <- bases[lag]
  lag.nms <- names (bases) <- rownames(mtf)[lag]
  nlag <- colSums(mtf[lag, , drop = FALSE])
  if (any(nlag > 1L))
    stop ("Interactions between lag terms are not allowed")

  ##
  lag.formula <- if (any(nlag == 0))
      reformulate(attr(drop.terms(mt, which(nlag == 0)), "term.labels"))
    else formula
  Bt <- as(t(model.matrix(lag.formula, mf.lag)), "sparseMatrix")
  if (attr(mt, "intercept"))
    Bt <- Bt[-1L, , drop = FALSE]
  rm (mf.lag)

  ## use text processing on the variable names to decide which
  ## groups of terms belong to which lag set (or not) and thus
  ## should be modeled as an independent random effect
  grp <- parse.names(lag.nms, rownames(Bt))
  if (any(grp == 0))  stop ("Unable to parse some basis terms. Please report")
  for (j in unique(grp)) {
    pn <- .names.func(j)
    data[[pn]] <- factor(rep(1:sum(grp == j), length.out = nrow(data)))
    fake.formula <- paste0(fake.formula, " + (1 | ", pn, ")")
    names (attr(grp, "dictionary"))[j] <- pn
  }
  bi <- if (length(lag.nms) == 1L)  rep(1L, length(attr(grp, "dictionary")))
    else  apply(sapply(lag.nms, grepl, attr(grp, "dictionary"), fixed = TRUE),
                1L, which)

  structure(
    list(formula = formula,
         lme4.formula = as.formula(fake.formula, env = environment(formula)),
         model = data,
         Bt = Bt,
         bases = bases,
         bi = bi,
         lag.group = grp
         ),
    class = "parsed.dlm"
    )
}
## interpret.dlm


