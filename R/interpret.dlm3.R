
#' @title Interpret a DLM Formula
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
#'
#' @details
#' Uses R's \code{\link[stats]{model.matrix}} mechanisms to build
#' and parse "random effects" components of distributed lag terms.
#'
#' @return
#' A \code{list} object...
#'

interpret.dlm3 <- function(formula, data,   ## random,
  .names.func = function(n) paste("pseudoGroups", n, sep = "")
) {
  fake.formula <- Reduce(paste, deparse(formula))
  mt <- attr(data, "terms")
  mtf <- attr(mt, "factor")
  ## attr(mt.tmp, "intercept") <- 0  ## RE design does not have intercept
  if (is.empty.model(mt))
    stop ("Empty model")

  ## to use R's own model.matrix() mechanisms for building main
  ## effects and interactions with DL terms, duplicate the model-frame
  ## but replace the "main effects" components of smooth DL terms with
  ## their "random effects" basis sets
  mf.lag <- data
  lag <- logical(ncol(mtf))  ## which terms in data are lag terms
  bases <- list()
  for (j in 1:ncol(mtf))  if (lag[j] <- is.SmoothLag(data[[j]])) {
    mf.lag[[j]] <- mf.lag[[j]]@random
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
  lag.formula <- reformulate(attr(drop.terms(mt, which(nlag == 0)), "term.labels"))
  Bt <- as(t(model.matrix(lag.formula, mf.lag)), "sparseMatrix")
  if (attr(mt, "intercept"))
    Bt <- Bt[-1L, , drop = FALSE]
  rm (mf.lag)

  fenms <- colnames (model.matrix(lme4::nobars(formula), mf))
  p.fenms <- suppressWarnings(parse.names(lag.nms, fenms))
  index <- list(covariates = p.fenms[character(1L)],
    smooth = lapply(lag.nms, function(i) p.fenms[i])
    )

  ## use text processing on the variable names to decide which
  ## groups of terms belong to which lag set (or not) and thus
  ## should be modeled as an independent random effect
  grp <- parse.names(lag.nms, rownames(Bt))
  for (j in unique(grp)) {
    pn <- .names.func(j)
    data[[pn]] <- factor(rep(1:sum(grp == j), length.out = nrow(data)))
    fake.formula <- paste0(fake.formula, " + (1 | ", pn, ")")
    names (index$smooth)[j] <- pn
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
         index = index,
         bi = bi,
         lag.names = lag.nms
         ),
    class = "parsed.dlm"
    )
}
## interpret.dlm




## Find and create index for unique names that have components
## matching base. The main idea here is that we want
## cr(lag, X):factor(group)1 to have a different index
## than cr(lag, X):factor(group)2, and so on
##
parse.names <- function(base, names) {
  regex <- sprintf("\\b?:?\\Q%s\\E[0-9.]*:?\\b?", base)
  base.index <- apply(sapply(regex, grepl, names), 1, which)
  if (is.list(base.index)) {
    base.index <- vapply(base.index, function(x) ifelse(length(x), x, NA),
                         numeric(1))
    warning ("Not all names have matching bases")
  }
  other <- sapply(strsplit(names, regex[base.index]), paste, collapse = ":")
  new.nms <- paste(base[base.index], other, sep = ":")
  new.nms[other == character(1L)] <- base[base.index[other == character(1L)]]
  new.nms[is.na(base.index)] <- character(1L)
  u <- unique(new.nms)
  structure(
    match(new.nms, u, nomatch = 0),
    dictionary = u,
    class = "pnames"
    )
}
## parse.names


"[.pnames" <- function(x, i) {
  ans <- which(x == which(attr(x, "dictionary") == i))
  if (length(ans)) ans else NA_integer_
}


## new.nms <- paste(ifelse(!is.na(base.index), base[base.index], ""),
##                  other, sep = ":")
##
## new.nms[is.na(base.index)] <- other[is.na(base.index)]
