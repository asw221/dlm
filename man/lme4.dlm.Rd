% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlm.R
\name{lme4.dlm}
\alias{lme4.dlm}
\title{Fit distributed lag models using lme4}
\usage{
lme4.dlm(parsed, family = gaussian(), control = list(), REML = FALSE, ...)
}
\arguments{
\item{parsed}{an interpreted \code{dlm} formula object returned by
\code{\link{interpret.dlm}}}

\item{family}{a description of the error distribution and link function to be used
in the model. The default is \code{gaussian(link = "identity")}.
See \code{stats::\link[stats]{family}} for possible family
functions and details}

\item{control}{either a \code{list} object with arguments to be passed to the
\code{lme4::\link[lme4]{lmerControl}} sequence, or the output of
\code{[g]lmerConrol} directly}

\item{REML}{if \code{TRUE} and a linear \code{\link{dlm}} model is specified,
\code{lme4.dlm} will use REML to fit the model. MLE will be
use otherwise}

\item{...}{other parameters to be passed to the \pkg{lme4}
\code{\link[lme4]{modular}} family functions}
}
\value{
an object that inherits from \code{lme4::\link[lme4]{merMod}}
  containing a fitted model
}
\description{
Fits an interpreted distributed lag model using \pkg{lme4}
\code{\link[lme4]{modular}} functions
}
\details{
Together with \code{\link{interpret.dlm}}, this function does the
main grunt work for \code{\link{dlm}}. Given an interpreted model,
\code{lme4.dlm} organizes the parsed data into the \pkg{lme4}
\code{\link[lme4]{modular}} functions to fit the model and returns
the fit as an \pkg{lme4} object.
}
