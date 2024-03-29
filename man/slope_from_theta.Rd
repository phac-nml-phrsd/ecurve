% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope_theta.R
\name{slope_from_theta}
\alias{slope_from_theta}
\title{Calculate the slope using theta.}
\usage{
slope_from_theta(theta)

slope_from_theta(theta)
}
\arguments{
\item{theta}{}
}
\value{
The slope parameter.

Numeric. The slope.
}
\description{
This is an internal function to conveniently parametrize
the slope such that the implied efficiency is naturally
constrainted between 0 and 1.
\code{Efficiency = 1 / (1 + exp(-theta))} and
\code{slope = 1 / log10(1+Efficiency)}
}
