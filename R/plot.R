
## plot.dlMod
## -------------------------------------------------------------------
#' @title Plot smoothed lag terms
#'
#' @description
#' Plot estimated lag coefficients
#'
#' @param x
#'   a fitted model object that inherits from \code{\link{dlMod}}
#' @param geom
#'   a \pkg{ggplot2} graph geometry. See Details
#' @param level
#'   a desired confidence interval level
#' @param scaled
#'   whether or not lag coefficients should be scaled
#'   (see \code{\link{estimands}})
#' @param ...
#'   additional graphical parameters to be passed to
#'   \code{ggplot2::\link[ggplot2]{facet_wrap}}
#'
#' @usage
#' ## S3 method for class 'dlMod'
#' plot(x, geom = c("pointrange", "line"), level = 0.95, scaled = TRUE, ...)
#'
#' @details
#' Plots estimated lag coefficients and confidence intervals to allow
#' convenient visual inspection of effects over different radii. Point
#' estimates and confidence intervals are computed using
#' \code{\link[=estimands]{confint}}, and plots are rendered using
#' minimal \pkg{ggplot2} commands so that the resulting graphics objects
#' are largely customizable.
#'
#' For now, the only supported plot geometries correspond to
#' \code{ggplot2::\link[ggplot2]{geom_pointrange}} and
#' \code{ggplot2::\link[ggplot2]{geom_line}}, and show the estimated
#' functions of Lag radii as either point and interval estimates,
#' or continuous functions with confidence bands, respectively. The
#' default option is \code{geom = "pointrange"} for point and interval
#' estimates.
#'
#' @return a \pkg{ggplot2} graphic object
#'
#' @name plotDlm
plotDlm <- function(x, y, level = 0.95, scaled = TRUE,
                    geom = c("pointrange", "line"), ...
) {
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

plot.dlMod <- plotDlm



## qqnorm.dlMod
## -------------------------------------------------------------------
qqnorm.dlMod <- function(y, ...) {
  r <- residuals(y, "pearson", scaled = TRUE)
  df <- data.frame(x = qnorm(ppoints(length(r))), y = sort(r))
  ggplot2::ggplot(df, aes(x = x)) +
    ggplot2::geom_point(aes(y = y), ...) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "darkgray") +
    ggplot2::labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
         title = "Normal Q-Q Plot"
         )
}
## qqnorm.dlMod



## .lag.ggplot (not exported)
## -------------------------------------------------------------------
.lag.ggplot <- function(.data, .geom, ...) {
  g <- ggplot2::ggplot(.data, aes(x = Lag)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~ term, ...)
  if (.geom == "pointrange")
    g + ggplot2::geom_pointrange(aes(y = Estimates, ymin = low, ymax = high),
                        size = rel(0.2))
  else if (.geom == "line")
    g + ggplot2::geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.35) +
      ggplot2::geom_line(aes(y = Estimates))
}
## .lag.ggplot


