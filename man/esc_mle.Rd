% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esc_mle.R
\name{esc_mle}
\alias{esc_mle}
\title{Fit ESC Model Using MLE}
\usage{
esc_mle(esc_data, CI = 0.95)
}
\arguments{
\item{esc_data}{Data frame containing data to be used to fit the model. Must
contain a column named "concentrations" with known sample concentrations, and
a column named "cqs" with corresponding Cq values. Non-detects should be
encoded by a Cq value of NaN}

\item{CI}{Numeric. Width of the confidence interval.}
}
\value{
esc object representing fitted model
}
\description{
Fit ESC Model Using MLE
}
