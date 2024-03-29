% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conc_est.R
\name{conc_mle}
\alias{conc_mle}
\title{Estimate concentration from Cq values}
\usage{
conc_mle(cqs, model, approximate = TRUE)
}
\arguments{
\item{cqs}{Numeric vector of Cq values from sample replicates, non-detects
coded as \code{NaN}.}

\item{model}{\code{esc} object representing a fitted model to use for estimation.}

\item{approximate}{Logical. If \code{TRUE} (the default), a faster but potentially
less accurate approximation for the likelihood function will be used at high
concentrations.}
}
\value{
The maximum likelihood estimation of the concentrations.
}
\description{
Given a list of Cq values from a set of technical replicates and
a fitted ESC model, generates a maximum likelihood estimate
of the concentrations.
}
\examples{

esc_data = data.frame(
     concentrations = c(1,1,10, 10, 100, 500, 500),
     cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
 mod = esc_mle(esc_data)

 new.cqs = c( 35, NaN, 36)
 x = conc_mle(cqs = new.cqs, model = mod)
 print(x)

}
