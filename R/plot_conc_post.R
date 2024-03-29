#' Plot Credible Intervals Concentration
#'
#' Given a \code{conc_int} object, plots either the posterior
#' probability density function or the cumulative density function
#' for the concentration,
#' along with the MLE estimate and the credible interval encoded
#' in the \code{conc_int} object
#'
#' @param interval A \code{conc_int} object as returned by
#' the function \code{conc_interval()} containing
#' the concentration posterior distribution and credible interval.
#' @param type Type of distribution function to plot, either \code{"pdf"} or \code{"cdf"}.
#' @param title String. Title for the plot. If unspecified, defaults to "Concentration
#' Posterior PDF" or "Concentration Posterior CDF", depending on the value of \code{type}
#'
#' @return A \code{ggplot} object representing generated plot.
#'
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' esc_data = data.frame(
#' concentrations = c(1,1,10, 10, 100, 500, 500),
#' cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
#' mod = esc_mle(esc_data)
#' new.cqs = c( 35, NaN, 36)
#' x = conc_interval(cqs = new.cqs, model = mod)
#' print(x$interval)
#' plot_conc_post(interval = x, type = 'pdf')
#' plot_conc_post(interval = x, type = 'cdf')
#'
plot_conc_post <- function(interval, type, title =
                            paste("Concentration Posterior", toupper(type))) {

  # Input checks

  if(!inherits(interval, 'conc_int')) {stop("interval is not a conc_int object")}
  if(type != "pdf" & type != "cdf") {stop("invalid type ", type, " specified")}
  if(!is.character(title) | !length(title) == 1) {stop("title must be a string")}

  # Compute positions of interval markers

  xvals <- unlist(interval$interval)[1:3]
  grid <- interval$distribution$concentration
  coords <- findInterval(xvals, grid)
  fun <- if(type == "pdf") {
    interval$distribution$pdf}
  else {
    interval$distribution$cdf}

  yvals <- (fun[coords] * (grid[coords + 1] - xvals) +
              fun[coords + 1] * (xvals - grid[coords])) / (grid[coords + 1] -
                                                             grid[coords])

  mle_label <- paste0("MLE: ", signif(interval$interval$mle, 3), "gc/rxn")

  interval_label <- paste0(round(interval$interval$level * 100), "% CI: ",
                           signif(interval$interval$lower, 3), " - ",
                           signif(interval$interval$upper, 3), "gc/rxn")

  point_df <- data.frame(x = xvals,
                         y = yvals,
                         Legend = c(interval_label,
                                    mle_label,
                                    interval_label))

  # Generate plot

  if(type == "pdf") {
    fun_mapping <- ggplot2::aes(x = .data$concentration, y = .data$pdf)
    label <- ggplot2::labs(x = "Concentration (gc/rxn)",
                           y = "Probability Density")
  }
  else {
    fun_mapping <- ggplot2::aes(x = .data$concentration, y = .data$cdf)
    label <- ggplot2::labs(x = "Concentration (gc/rxn)",
                           y = "Cumulative Probability")
  }
  plot <- ggplot2::ggplot(interval$distribution, fun_mapping) +
    label +
    ggplot2::geom_line() +
    ggplot2::geom_segment(
      data = point_df,
      mapping = ggplot2::aes(x     = .data$x,
                             y     = 0,
                             xend  = .data$x,
                             yend  = .data$y,
                             color = .data$Legend))
  return(plot)
}
