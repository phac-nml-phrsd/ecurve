
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
#' constrainted between 0 and 1.
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

mylogistic <- function(x, minvalue, maxvalue) {
  return( minvalue + maxvalue / (1 + exp(-x)))
}

mylogistic_inverse <- function(y, minvalue, maxvalue) {
  return( -log( maxvalue / (y - minvalue) -1 )  )
}

# xx = -20:20
# yy = sapply(xx, mylogistic, minvalue = 2, maxvalue = 35)
# zz = sapply(yy, mylogistic_inverse, minvalue = 2, maxvalue = 35)
#
# check = mylogistic_inverse(mylogistic(x=xx, minvalue = 2, maxvalue = 35), minvalue = 2, maxvalue = 35)
# plot(xx,check);abline(a=0,b=1)
# plot(xx,yy, typ='o')
# plot(yy,zz, typ='o')
