

## Non-Exported Misc Utility Functions for builtenvir Package
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





## Find and create index for unique names that have components
## matching base. The main idea here is that we want
## cr(lag, X):factor(group)1 to have a different index
## than cr(lag, X):factor(group)2, and so on
##
parse.names <- function(base, names, .warn = TRUE) {
  regex <- sprintf("\\b?:?\\Q%s\\E[0-9.]*:?\\b?", base)
  base.index <- apply(sapply(regex, grepl, names), 1, which)
  if (is.list(base.index)) {
    base.index <- vapply(base.index, function(x) ifelse(length(x), x, NA),
                         numeric(1))
    if (.warn)  warning ("Not all names have matching bases")
  }
  other <- sapply(strsplit(names, regex[base.index]), paste, collapse = ":")
  new.nms <- paste(base[base.index], other, sep = ":")
  new.nms[other == character(1L)] <- base[base.index[other == character(1L)]]
  new.nms[is.na(base.index)] <- character(1L)
  u <- unique(new.nms)
  structure(
    match(new.nms, u, nomatch = 0),
    Dictionary = u,
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


## new.nms <- paste(ifelse(!is.na(base.index), base[base.index], ""),
##                  other, sep = ":")
##
## new.nms[is.na(base.index)] <- other[is.na(base.index)]

