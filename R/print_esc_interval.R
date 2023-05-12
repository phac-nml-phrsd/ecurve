#' @title Print ESC interval.
#'
#' @param x A ESC interval object as returned by the function \code{conc_interval()}.
#'
#' @return Print out of the maximum likelihood estimator,
#' upper and lower bound of the pre-specified statistical confidence level.
#'
#' @export
#'
#' @examples
#'
#' esc_data = data.frame(
#'    concentrations = c(1,1,10, 10, 100, 500, 500),
#'    cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
#' model = esc_mle(esc_data)
#' cqs = new.cqs = c( 34, NaN, 36)
#' x = conc_interval(cqs = new.cqs, model = model)
#' print_esc_interval(x)
#'
#'
print_esc_interval <- function(x) {

  if(!inherits(x, 'conc_int')){
    cat('The object is not a ESC `conc_int` object, print_esc_interval()` cannot print it.')
    warning('Not a `conc_int` object : nothing to print.')
    return()
  }

  lvl = x$interval$level
  lo  = x$interval$lower
  mle = x$interval$mle
  hi  = x$interval$upper

  cat(paste0(
    '---- ESC interval ----\n\n',
    'Confidence level = ', lvl, '\n',
    'upper bound      = ', hi, '\n',
    'max. likelihood  = ', mle, '\n',
    'lower bound      = ', lo, '\n',
    '----------------------\n'
  ))

}
