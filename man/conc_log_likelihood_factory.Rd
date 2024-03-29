% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conc_est.R
\name{conc_log_likelihood_factory}
\alias{conc_log_likelihood_factory}
\title{Factory Function to Create Concentration Log Likelihood Function}
\usage{
conc_log_likelihood_factory(cqs, model)
}
\arguments{
\item{cqs}{Numeric vector of observed Cq values, with non-detects coded as NaN}

\item{model}{\code{esc} object representing fitted model to use for estimation.}
}
\value{
Function that takes single input concentration and returns negative of
corresponding log likelihood
}
\description{
Given ESC Model parameters and observed Cq data, returns function implementing
negative log likelihood calculation for concentration estimation
}
