

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


## Return arbitrary dimension of array if it exists; otherwise return 1
NDIM <- function(x, n = 1L)
  if (length(d <- dim(x)) > n - 1L) d[n] else 1L




