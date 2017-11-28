
## Methods for SmoothLag Class
## -------------------------------------------------------------------

## getGeneric("[")
## ## standardGeneric for "[" defined from package "base"
## ## function (x, i, j, ..., drop = TRUE)
## ## standardGeneric("[", .Primitive("["))

## getGeneric("[<-")
## ## standardGeneric for "[<-" defined from package "base"
## ## function (x, i, j, ..., value)
## ## standardGeneric("[<-", .Primitive("[<-"))

## x[]
setMethod("[", signature = c(x = "SmoothLag", i = "missing", j = "missing",
                             drop = "ANY"),
          function(x, i, j, ..., drop) x
          )


setMethod("[", signature = c(x = "SmoothLag", i = "numeric", j = "missing",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            SmoothLag(x@.Data[i, , drop = FALSE],
                      basis = x@basis,
                      random = x@random[i, , drop = FALSE],
                      signature = x@signature
                      )
          })




setMethod("[", signature = c(x = "SmoothLag", i = "missing", j = "numeric",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            x@.Data <- x@.Data[, j, drop = FALSE]
            x
          })


setMethod("[", signature = c(x = "SmoothLag", i = "numeric", j = "numeric",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            SmoothLag(x@.Data[i, j, drop = FALSE],
                      basis = x@basis,
                      random = x@random[i, , drop = FALSE],
                      signature = x@signature
                      )
          })


setMethod("[", signature = c(x = "SmoothLag", i = "logical", j = "missing",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            SmoothLag(x@.Data[i, , drop = FALSE],
                      basis = x@basis,
                      random = x@random[i, , drop = FALSE],
                      signature = x@signature
                      )
          })

setMethod("[", signature = c(x = "SmoothLag", i = "missing", j = "logical",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            x@.Data <- x@.Data[, j, drop = FALSE]
            x
          })

setMethod("[", signature = c(x = "SmoothLag", i = "logical", j = "logical",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            SmoothLag(x@.Data[i, j, ..., drop = FALSE],
                      basis = x@basis,
                      random = x@random[i, , drop = FALSE],
                      signature = x@signature
                      )
          })

setMethod("[", signature = c(x = "SmoothLag", i = "numeric", j = "logical",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            SmoothLag(x@.Data[i, j, ..., drop = FALSE],
                      basis = x@basis,
                      random = x@random[i, , drop = FALSE],
                      signature = x@signature
                      )
          })

setMethod("[", signature = c(x = "SmoothLag", i = "logical", j = "numeric",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            SmoothLag(x@.Data[i, j, ..., drop = FALSE],
                      basis = x@basis,
                      random = x@random[i, , drop = FALSE],
                      signature = x@signature
                      )
          })


## alternative catch-all
setMethod("[", signature = c(x = "SmoothLag", i = "ANY", j = "ANY",
                             drop = "ANY"),
          function(x, i, j, ..., drop)
            stop ("invalid or not yet implemented dimension set")
          )
