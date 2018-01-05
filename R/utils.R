

## Exported misc utility functions
## -------------------------------------------------------------------

## Find and create index for unique names that have components
## matching base. The main idea here is that we want
## cr(lag, X):factor(group)1 to have a different index
## than cr(lag, X):factor(group)2, and so on
##
#' @title Partition DL term names
#'
#' @description
#' Parse variable names and create a unique integer label for each unique use
#' of a matching base case. \code{parse.names} is used internally by several
#' \pkg{dlmBE} functions to determine which distributed lag terms
#' share the same random effects variance/penalty.
#'
#' @param base
#'   a character vector with literal matches in \code{names}
#' @param names
#'   a character vector with substrings that are literal matches to
#'   elements in \code{base}
#'
#' @return
#' an S3 object of class \code{"pnames"} which is of \code{length}
#' equal to \code{length(names)} and stores a unique integer label
#' for each matching case in \code{base} and a look up dictionary
#' for each integer.
#'
#' Indexing in \code{pnames} objects takes a character string argument
#' and returns the vector indices of the matching \code{base} uses.
#' See Examples.
#'
#' @examples
#' lg <- 1:5
#' Z <- matrix(rbinom(60, 10, 0.1), 12)
#' group <- rep(1:4, each = 3)
#' base <- "cr(lg, Z)"
#' names <- colnames(model.matrix(~ cr(lg, Z):factor(group)))
#' pn <- parse.names(base, names, FALSE)
#' names
#' pn
#' pn[""]
#' pn["cr(lg, Z):factor(group)2"]
#'
parse.names <- function(base, names, .warn = TRUE) {
  regex <- sprintf(":?\\b?\\Q%s\\E[0-9.]*\\b?:?", base)
  base.index <- apply(sapply(regex, grepl, names), 1, which)
  if (is.list(base.index)) {
    base.index <- vapply(base.index, function(x) ifelse(length(x), x, NA),
                         numeric(1))
    if (.warn)  warning ("Not all names have matching bases")
  }
  other <- sapply(strsplit(names, regex[base.index]),
                  function(x) paste(x[nzchar(x)], collapse = ":"))
  new.nms <- paste(base[base.index], other, sep = ":")
  new.nms[!nzchar(other)] <- base[base.index[!nzchar(other)]]
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
  .find <- function(x, i) which(x == which(attr(x, "dictionary") == i))
  ans <- unlist(lapply(i, function(j) .find(x, j)))
  if (!is.null(ans)) ans else NA_integer_
}
## [.pnames







## Non-Exported misc utility functions
## -------------------------------------------------------------------

## Throw error for argument of unrecognized type
.Unrecognized <- function(name, class, .n = 1) {
  caller.name <- as.character(sys.call(.n))[1]
  stop (caller.name, ": Unrecognized type for '",
        name, "': ", paste(class, collapse = ", "),
        call. = FALSE
        )
}
## .Unrecognized


## Warn if arguments dropped
.Ignored <- function(..., call. = FALSE, immediate. = FALSE,
                     noBreaks. = FALSE, domain = NULL
                     ) {
  l... <- list(...)
  if (length(l...))
    warning ('Arguments, "', paste(names(l...), collapse = ", "), '" ignored',
             call. = call., immediate. = immediate., noBreaks. = noBreaks.,
             domain = domain
             )
}





