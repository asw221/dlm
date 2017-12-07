


## plot.dlMod
## -------------------------------------------------------------------
#' @title Plot smoothed lag terms
#'
#' @description
#' Plot estimated lag-spline effects
#'
#' @usage
#' plot(x, geom = c("pointrange", "line"), level = 0.95, scaled = TRUE, ...)
#'
#' @return a \code{ggplot2} graphic object
#'
#' @name plotDlm
plot.dlMod <- function(x, y, geom = c("pointrange", "line"),
                       level = 0.95, scaled = TRUE, ...) {
  if (!missing(y) && missing(geom)) geom <- y
  else if (!missing(y)) .Ignored(y = y)
  geom <- match.arg(geom)
  ci <- confint(x, level = level, scaled = scaled)
  ci <- ci[, c(1:2, ncol(ci))]
  colnames (ci) <- c("Estimates", "low", "high")
  li <- lagIndex(x)
  df <- do.call("rbind", lapply(li,
          function(i) data.frame(ci[i, ], check.names = FALSE))
          )
  df[["Lag"]] <- unlist(lapply(x@index, function(i) x@bases[[i]]@x))
  df[["term"]] <- rep(names(li), lengths(li))
  .lag.ggplot(df, geom, ...)
}
## plot.dlMod



## qqnorm.dlMod
## -------------------------------------------------------------------
qqnorm.dlMod <- function(y, ...) {
  r <- residuals(y, "pearson", scaled = TRUE)
  df <- data.frame(x = qnorm(ppoints(length(r))), y = sort(r))
  ggplot(df, aes(x = x)) +
    geom_point(aes(y = y), ...) +
    geom_abline(intercept = 0, slope = 1, color = "darkgray") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
         title = "Normal Q-Q Plot"
         )
}
## qqnorm.dlMod



## .lag.ggplot (not exported)
## -------------------------------------------------------------------
.lag.ggplot <- function(.data, .geom, ...) {
  g <- ggplot(.data, aes(x = Lag)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ term, ...)
  if (.geom == "pointrange")
    g + geom_pointrange(aes(y = Estimates, ymin = low, ymax = high),
                        size = rel(0.2))
  else if (.geom == "line")
    g + geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.35) +
      geom_line(aes(y = Estimates))
}
## .lag.ggplot


