



plot.dlMod <- function(x, y, geom = c("pointrange", "line"),
                       level = 0.95, ...) {
  geom <- match.arg(geom)
  ci <- confint(x, level = level, ...)
  ci <- ci[, c(1, ncol(ci))]
  li <- lagIndex(x)
  df <- do.call("rbind", lapply(li,
          function(i) data.frame(ci[i, ], check.names = FALSE))
          )
  df[["term"]] <- rep(names(li), lengths(li))

  invisible (0)
}

