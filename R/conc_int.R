#' Generate New conc_int Object
#'
#' Creates new object of class conc_int, containing a list specifying the upper
#' and lower endpoints of the credible interval,the MLE concentration
#' estimate, and the level of the interval, as well as a data frame containing
#' the PDF and CDF of posterior concentration distribution evaluated at a grid
#' of points
#'
#' @param mle MLE estimate of concentration
#' @param interval Numeric vector representing credible interval: first element
#' is the lower endpoint, second is the upper endpoint
#' @param grid Numeric vector containing grid of value at with PDF/CDF were
#' evaluated
#' @param pdf Numeric vector of PDF values evaluated at each point in grid
#' @param cdf Numeric vector of CDF values evaluated at each point in grid
#' @param level Numeric. Level of the credible interval
#'
#' @return Object of class conc_int
#'
#' @examples
#'
new_conc_int <- function(mle, interval, grid, pdf, cdf, level) {
  norm_factor <- cdf[length(cdf)]
  structure(list(
    interval = list(lower = interval[1], mle = mle, upper = interval[2], level = level),
    distribution = data.frame(concentration = grid, pdf = pdf/norm_factor,
                              cdf = cdf/norm_factor)
  ), class = "conc_int")
}
