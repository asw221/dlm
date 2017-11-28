
#' Simulated Built Environment Data
#'
#' @format
#' The first column (\code{Y}) is the outcome variable, covariates \code{Gender}
#' and \code{Age} come next, and are followed by 100 location description
#' variables. Each column of these 100 location count variables corresponds to
#' a distance lag equal to the values in \code{seq(0.1, 10, 0.1)}.
#'
#' @docType data
#'
#' @usage data(simdata)
#'
#' @examples
#' data(simdata)
#' simdata[1:10, 1:5]
#'
#' @inherit dlm references
#'
"simdata"
