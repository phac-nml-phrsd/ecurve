
#' @title Calculate the slope using theta.
#'
#' @description This is an internal function to conveniently parametrize
#' the slope such that the implied efficiency is naturally
#' constrainted between 0 and 1.
#' `Efficiency = 1 / (1 + exp(-theta))` and
#'  `slope = 1 / log10(1+Efficiency)`
#'
#' @param theta Numeric.
#'
#' @return The slope parameter.
#'
slope_from_theta <- function(theta) {
  tmp = 1 + 1 / (1 + exp(-theta))
  return( -1 / log10(tmp) )
}


#' @title Calculate theta using the slope .
#'
#' @description This is an internal function to conveniently parametrize
#' the slope such that the implied efficiency is naturally
#' constrained between 0 and 1.
#' `Efficiency = 1 / (1 + exp(-theta))` and
#'  `slope = 1 / log10(1+Efficiency)`
#'
#' @param slope Numeric.
#'
#' @return The theta parameter.
#'
theta_from_slope <- function(slope) {
  tmp  = 10^(1/-slope)
  tmp1 = tmp - 1
  tmp2 = 2 - tmp
  return(log(tmp1/tmp2))
}

#' Logistic function.
#'
#' @param x Numeric
#' @param minvalue Minimum value.
#' @param maxvalue Maximum value.
#'
#' @return logistic(x)
#'
mylogistic <- function(x, minvalue, maxvalue) {
  return( minvalue + maxvalue / (1 + exp(-x)))
}

#' Inverse of the logistic function.
#'
#' @param y Numeric
#' @param minvalue Minimum value.
#' @param maxvalue Maximum value.
#'
#' @return logistic^-1(y)
#'
mylogistic_inverse <- function(y, minvalue, maxvalue) {
  return( -log( maxvalue / (y - minvalue) -1 )  )
}
