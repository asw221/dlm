
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

._interpret.dlm <- function(formula, data,   ## random,
  .names.func = function(n) paste("pseudoGroups", n, sep = "")
) {
  fake.formula <- Reduce(paste, deparse(formula))
  mt <- mt.tmp <- attr(data, "terms")
  ## attr(mt.tmp, "intercept") <- 0  ## RE design does not have intercept
  if (is.empty.model(mt))
    stop ("Empty model")

  ## to use R's own model.matrix() mechanisms for building main
  ## effects and interactions with DL terms, duplicate the model-frame
  ## but replace the "main effects" components of smooth DL terms with
  ## their "random effects" basis sets
  mf.lag <- data
  lag <- logical(ncol(data))  ## which terms in data are lag terms
  bases <- list()
  for (j in 1:ncol(data))  if (lag[j] <- inherits(data[[j]], "SmoothLag")) {
    mf.lag[[j]] <- mf.lag[[j]]@random
    bases[[j]] <- data[[j]]@basis
    data[[j]] <- data[[j]]@.Data
  }
  if (!any(lag))
    stop ("Model has no apparent DL terms")

  bases <- bases[lag]
  mtf <- attr(mt, "factor")
  lag.nms <- names(bases) <- rownames(mtf)[lag]
  nlag <- colSums(mtf[lag, , drop = FALSE])
  if (any(nlag > 1))
    stop ("Interactions between lag terms are not allowed")

  is.lag <- as.logical(nlag)
  attr(mt.tmp, "factors") <- mtf[, is.lag, drop = FALSE]
  attr(mt.tmp, "term.labels") <- attr(mt.tmp, "term.labels")[is.lag]
  attr(mt.tmp, "order") <- attr(mt.tmp, "order")[is.lag]

  ## use text processing on the variable names to decide which
  ## groups of terms belong to which lag set (or not) and thus
  ## should be modeled as an independent random effect
  X <- model.matrix(mt, data)
  Z <- model.matrix(mt.tmp, mf.lag)[, -1]
  P <- ncol(X)
  ## nms <- c(colnames(Z), colnames(X))
  nms <- c(colnames(X), colnames(Z))
  rm (X)

  grp <- parse.names(lag.nms, nms)  ## 0 for non-lag covariates
  index <- list(covariates = NULL, smooth = list())

  ## Build "pseudo-groups"
  for (j in unique(grp)) {
    if (j == 0)
      index$covariates <- which(grp == j)
    else {
      index$smooth[[j]] <- which(grp == j)
      k <- sum(index$smooth[[j]] > P)
      pn <- .names.func(j)
      data[[pn]] <- factor(rep(1:k, length.out = nrow(data)))
      fake.formula <- paste0(fake.formula, " + (1 | ", pn, ")")
    }
  }
  names (index$smooth) <- attr(grp, "dictionary")
  bi <- if (length(lag.nms) == 1)  rep(1, length(attr(grp, "dictionary")))
    else  apply(sapply(lag.nms, grepl, attr(grp, "dictionary"), fixed = TRUE),
                1, which)

  ## ## add other random effects
  ## if (!missing(random) && length(random))
  ##   lme4.formula <- paste(lme4.formula,
  ##     paste(unlist(lapply(random, function(x) gsub("~", "", deparse(x), fixed = TRUE))),
  ##           collapse = " + "), sep = " + ")

  res <- list(formula = formula,
              lme4.formula = as.formula(fake.formula,
                                        env = environment(formula)),
              terms = mt,
              model = data,
              B = Z,
              np = c(nrow(data), P),
              index = index,
              bases = bases,
              bi = bi
              )
  class (res) <- "parsed.dlm"
  return (res)
}
## interpret.dlm




## Find and create index for unique names that have components
## matching base. The main idea here is that we want
## cr(lag, X):factor(group)1 to have a different index
## than cr(lag, X):factor(group)2, and so on
parse.names <- function(base, names) {
  regex <- sprintf("%s[^:]*", gsub("([.\\\\|()[{^$*+?])", "\\\\\\1", base))
  u <- character(0)
  for (i in 1:length(base)) {
    ndx <- grep(base[i], names, fixed = TRUE)
    names[ndx] <- gsub(regex[i], base[i], names[ndx])
    u <- c(u, unique(names[ndx]))
  }
  m <- match(names, u, nomatch = 0)
  attr(m, "dictionary") <- u
  return (m)
}





## ## deviance:
## if (!is.na(dev <- fit@devcomp$cmp["REML"])) dev else fit@devcomp$cmp["dev"]
## ## either of these will be NA o/w
