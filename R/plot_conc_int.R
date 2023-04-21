#' Plot Credible Intervals Concentration
#'
#' Given a \code{conc_int} object, plots either the posterior
#' probability density function or the cumulative density function
#' for the concentration,
#' along with the MLE estimate and the credible interval encoded
#' in the \code{conc_int} object
#'
#' @param interval a \code{conc_int} object containing the concentration posterior distribution
#' and credible interval.
#' @param type distribution function to plot, either \code{"pdf"} or \code{"cdf"}.
#'
#' @return A \code{ggplot} object representing generated plot.
#'
#' @export
#'
#' @examples
#'
#'
plot_conc_int <- function(interval, type) {

  # Input checks
  if(!inherits(interval, 'conc_int')) {stop("interval is not a conc_int object")}
  if(type != "pdf" & type != "cdf") {stop("invalid type ", type, " specified")}

  #Compute positions of interval markers
  xvals <- unlist(interval$interval)
  grid <- interval$distribution$concentration
  coords <- findInterval(xvals, grid)
  fun <- if(type == "pdf") {interval$distribution$pdf} else {interval$distribution$cdf}
  yvals <- (fun[coords] * (grid[coords + 1] - xvals) +
              fun[coords + 1] * (xvals - grid[coords])) / (grid[coords + 1] -
                                                             grid[coords])
  point_df <- data.frame(x = xvals,
                         y = yvals,
                         Legend = c("Interval", "MLE", "Interval"))

  # Generate plot
  if(type == "pdf") {
    fun_mapping <- ggplot2::aes(x = concentration, y = pdf)
    label <- ggplot2::ylab("Probability Density")
  }
  else {
    fun_mapping <- ggplot2::aes(x = concentration, y = cdf)
    label <- ggplot2::ylab("Cumulative Probability")
  }
  plot <- ggplot2::ggplot(interval$distribution, fun_mapping) + label +
    ggplot2::geom_line() +
    ggplot2::geom_segment(data = point_df,
                          mapping = ggplot2::aes(x = x, y = 0, xend = x,
                                                 yend = y, color = Legend))
 return(plot)
}
