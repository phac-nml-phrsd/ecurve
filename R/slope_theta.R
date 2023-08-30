
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


#' Calculate the regression slope from theta
#'
#' @param theta
#'
#' @return Numeric. The slope.

slope_from_theta <- function(theta) {
  eff = 1 / (1 + exp(-theta))
  slope = -1 / log10(1 + eff)
  return(slope)
}

#' Logistic function.
#'
#' @param x Numeric
#' @param minvalue Minimum value.
#' @param maxvalue Maximum value.
#' @param shallowness Numeric. Determines how steep the logistic function is.
#'
#' @return logistic(x)
#'
mylogistic <- function(x, minvalue, maxvalue, shallowness = 25) {
  return( minvalue + (maxvalue - minvalue) / (1 + exp(-x/shallowness)))
}

#' Inverse of the logistic function.
#'
#' @param y Numeric
#' @param minvalue Minimum value.
#' @param maxvalue Maximum value.
#' @param shallowness Numeric. Determines how steep the logistic function is.
#'
#' @return logistic^-1(y)
#'
mylogistic_inverse <- function(y, minvalue, maxvalue, shallowness = 25) {
  return( -log( (maxvalue - minvalue) / (y - minvalue) -1 ) * shallowness  )
}


if(0){

  minvalue = 10
  maxvalue = 60
  x = seq(minvalue - 5,maxvalue+5,by=0.2)
  y = mylogistic_inverse(x,minvalue,maxvalue)
  z = mylogistic(y,minvalue,maxvalue)
  par(mfrow=c(1,3))
  plot(x,z) ; abline(0,1, col='red')
  plot(x,y,typ='l') ; grid()
  plot(y,z,typ='l') ; grid()

}


