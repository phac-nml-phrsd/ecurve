% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_esc.R
\name{plot_esc_model}
\alias{plot_esc_model}
\title{Plot Fitted ESC Model}
\usage{
plot_esc_model(
  model,
  PI = 0.95,
  title = "ESC model fit",
  xlimits = NULL,
  approximate = TRUE
)
}
\arguments{
\item{model}{\code{esc} object representing fitted model.}

\item{PI}{Numeric. Width of the probability interval (must be between 0 and 1).}

\item{title}{String. The title of the plot. If unspecified, defaults to "ESC
model fit"}

\item{xlimits}{Optional numeric vector of length 2. If provided, passed to
\code{ggplot2} to specify x-axis limits of plot}

\item{approximate}{Logical. If \code{TRUE} (the default), a faster but potentially
less accurate approximation for the likelihood function will be used at high
concentrations.}
}
\value{
A \code{ggplot} object representing plot.
}
\description{
Plots the fitted median and probability interval of
an ESC model Cq vs concentration.
}
\examples{
esc_data = data.frame(
    concentrations = c(1,1,10, 10, 100, 500, 500),
    cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
mod = esc_mle(esc_data)
plot_esc_model(mod)


}
