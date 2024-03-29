% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_esc_interval.R
\name{print_esc_interval}
\alias{print_esc_interval}
\title{Print ESC interval.}
\usage{
print_esc_interval(x)
}
\arguments{
\item{x}{A ESC interval object as returned by the function \code{conc_interval()}.}
}
\value{
Print out of the maximum likelihood estimator,
upper and lower bound of the pre-specified statistical confidence level.
}
\description{
Print ESC interval.
}
\examples{

esc_data = data.frame(
   concentrations = c(1,1,10, 10, 100, 500, 500),
   cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
model = esc_mle(esc_data)
cqs = new.cqs = c( 34, NaN, 36)
x = conc_interval(cqs = new.cqs, model = model)
print_esc_interval(x)


}
