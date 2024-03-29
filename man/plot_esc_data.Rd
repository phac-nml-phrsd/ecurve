% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_esc.R
\name{plot_esc_data}
\alias{plot_esc_data}
\title{Plot Standard Curve Calibration Data}
\usage{
plot_esc_data(esc_data, xlimits = NULL)
}
\arguments{
\item{esc_data}{Data frame containing data to bt used to fit the model. Must
contain a column named \code{concentrations} with known sample concentrations,
and a column named \code{cqs} with corresponding Cq values.
Non-detects should be encoded by a Cq value of \code{NaN}.}

\item{xlimits}{Optional numeric vector of length 2. If provided, passed to
\code{ggplot2} to specify x-axis limits of plot}
}
\value{
A \code{ggplot} object representing plot.
}
\description{
Given a data frame containing concentrations
and corresponding Cq values, plots Cq values against
log-transformed concentrations.
}
\examples{
esc_data = data.frame(
    concentrations = c(1,1,10, 10, 100, 500, 500),
    cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
plot_esc_data(esc_data)

}
