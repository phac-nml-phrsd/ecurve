#' Plot Concentration Credible Intervals
#'
#' Given conc_int object, plots either posterior pdf or cdf for concentration,
#' along with MLE estimate and credible interval encoded in conc_int object
#'
#' @param interval conc_int object containing concentration posterior distribution
#' and credible interval
#' @param type distribution function to plot, either "pdf" or "cdf"
#'
#' @return ggplot object representing generated plot
#' @export
#'
#' @examples
plot_conc_int <- function(interval, type) {

  #Input checks
  if(class(interval) != "conc_int") {stop("interval is not a conc_int object")}
  if(type != "pdf" & type != "cdf") {stop("invalid type ", type, " specified")}

  #Compute positions of interval markers
  xvals <- unlist(interval$interval)
  grid <- interval$distribution$concentration
  coords <- findInterval(xvals, grid)
  fun <- if(type == "pdf") {interval$distribution$pdf} else {interval$distribution$cdf}
  yvals <- (fun[coords] * (grid[coords + 1] - xvals) +
              fun[coords + 1] * (xvals - grid[coords])) / (grid[coords + 1] -
                                                             grid[coords])
  point_df <- data.frame(x = xvals, y = yvals, Legend = c("Interval", "MLE", "Interval"))

  #Generate plot
  if(type == "pdf") {fun_mapping <- ggplot2::aes(x = concentration, y = pdf)}
  else {fun_mapping <- ggplot2::aes(x = concentration, y = cdf)}
  plot <- ggplot2::ggplot(interval$distribution, fun_mapping) +
    ggplot2::geom_line() +
    ggplot2::geom_segment(data = point_df,
                          mapping = ggplot2::aes(x = x, y = 0, xend = x,
                                                 yend = y, color = Legend))

  #Diplay and return plot
  print(plot)
  plot
}
