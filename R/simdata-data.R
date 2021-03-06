
#' Simulated built environment data
#'
#' @docType data
#' @usage data(simdata)
#'
#' @format
#' The first column (\code{Y}) is the outcome variable (intended to mimic
#' negative BMI score), covariates \code{Gender} (binary)
#' and \code{Age} (years) come next, and are followed by 100 location description
#' variables each corresponding to
#' a lag radius equal to the values in \code{seq(0.1, 10, 0.1)}. These
#' 100 location variables store the counts of simulated built environment
#' features within each radius for each participant.
#'
"simdata"
