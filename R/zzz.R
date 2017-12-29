
#' dlm: Distributed lag models in R using lme4
#'
#' @inherit dlm references
#' @inherit basis references
#'
#' @description
#' This package fits the distributed lag models (DLMs) described by Baek et
#' al (2016) and Baek et al (2017), which estimate the association between
#' the presence of built environment features and an outcome as a function
#' of distance between the locations for study participants and locations for
#' environment features or community resources. These models circumvent the
#' need to pre-specify a radius within which to measure the availability of
#' community resources.  Distributed lag models have a long history in a
#' variety of fields. For built environment research, we define the lagged
#' exposure as the value of an environment feature between two radii,
#' \eqn{r_{l-1}} and \eqn{r_l}
#' from study locations, \eqn{l = 1, 2, \ldots, L}{l = 1, 2, ..., L},
#' where \eqn{r_0 = 0}; e.g., the lagged exposure is the number
#' of convenience stores within “ring”-shaped areas around study participants
#' residential address.  The package supports generalized linear regression
#' models, as well as generalized linear mixed models.  In both instances,
#' multiple lagged exposure covariates maybe included, as well as
#' interactions between the lagged covariates and other categorical
#' covariates (e.g., quartiles of age).
#'
#' Let \eqn{Y_{ij}} be an outcome measured at location \eqn{i} at visit
#' \eqn{j}, and \eqn{X_{ij}(r_{l-1}; r_l)}
#' be an environment feature measured during visit \eqn{j} within a
#' ring-shaped area around location \eqn{i} between radii \eqn{r_{l-1}}
#' and \eqn{r_l}; and \eqn{r_L} be the maximum distance
#' around locations beyond which there is no association between the
#' environment feature and the outcome. A typical unadjusted generalized
#' linear mixed model that can be fitted in this version of the package is,
#'
#' \deqn{g(E(Y_{ij} \vert b_i)) = \beta_0 + \sum_{l=1}^L \beta(r_{l-1}; r_l) X(r_{l-1}; r_l) + W_{ij} b_i}{g(E(Y_{ij} | b_i)) = \beta_0 + \sum_{l=1}^L \beta(r_{l-1}; r_l) * X(r_{l-1}; r_l) + W_{ij} * b_i}
#'
#' where \eqn{g(\cdot)}{g()} is a link function appropriate for the
#' distribution of the outcome; \eqn{\beta_0}
#' represents an intercept; the association of the environment feature
#' measured between radii \eqn{r_{l-1}} and \eqn{r_l} and the outcome is
#' \eqn{\beta(r_{l-1}; r_l)}; and \eqn{W_{ij}} are covariates related
#' to random effects, \eqn{b_i} (e.g., random intercepts and slopes).
#' The coefficients \eqn{\beta(r_{l-1}; r_l)} are constrained to follow a
#' smooth function of distance
#' from the locations of interest; the constraint is imposed by modeling the
#' coefficients using smoothing splines. Other models could be used, although
#' smoothing splines are the only supported option at this time.
#'
#' The model easily simplifies to generalized linear regression modes
#' (e.g., when there is only one visit), and can be extended in the following
#' directions. Adjustment covariates can be easily included. In addition,
#' interaction terms between covariates and the DL covariates are
#' also supported. For example, terms such as:
#' \eqn{\sum_{l=1}^L \theta(r_{l-1}; r_l) X(r_{l-1}; r_l) Z_i}{\sum_{l=1}^L \theta(r_{l-1}; r_l) * X(r_{l-1}; r_l) * Z_i},
#' where \eqn{Z_i} is another covariate, can be included. The interaction
#' coefficients \eqn{\theta(r_{l-1}; r_l)}  have the
#' usual interpretation, but the magnitude of the interaction can vary over
#' distance from locations of interest; \eqn{\theta(r_{l-1}; r_l)} are also
#' constrained using smoothing splines. Finally, weighted regression models
#' are also supported.
#'
#' We assume the user has calculated distances from every participant’s
#' location to every community resource/feature. The distances can be network
#' distances or Euclidian distances. Those distances are then used to calculate
#' the distributed lag covariates, \eqn{X(r_{l-1}; r_l)}, by specifying \eqn{L}
#' and the radii \eqn{r_l}, \eqn{l = 1, 2, \ldots, L}{l = 1, 2, ..., L}.
#' See Baek et al (2016) for guidance on choosing \eqn{L} and \eqn{r}.
#'
#' The package includes a series of functions to pass formulas and data to
#' \pkg{lme4}, which is used for estimation of the DLM. All those functions are
#' documented in this manual, although a typical user will primarily interact
#' with XXX, xxx, and xxx. For example:
#'
#'
#' @docType package
#' @name dlm-package
NULL
