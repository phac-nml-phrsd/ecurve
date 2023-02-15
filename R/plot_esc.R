#' Plot Standard Curve Calibration Data
#'
#' Given data frame containing concentrations and corresponding Cq values, plots
#' Cq values against log-transformed concentrations
#'
#' @param esc_data Data frame containing data to bt used to fit the model. Must
#' contain a column named "concentrations" with known sample concentrations, and
#' a column named "cqs" with corresponding Cq values. Non-detects should be
#' encoded by a cq value of NaN
#'
#' @return ggplot object representing plot
#' @export
#'
#' @examples
plot_esc_data <- function(esc_data) {

  # --- Input checks

  if(!is.data.frame(esc_data)) {stop("esc_data must be a data frame")}
  if(!"concentrations" %in% names(esc_data)) {
    stop("esc_data must contain concentrations column")
  }
  if(!"cqs" %in% names(esc_data)) {
    stop("esc_data must contain cqs column")
  }

  concentrations <- esc_data[,"concentrations"]
  cqs <- esc_data[,"cqs"]

  if(!all(is.numeric(concentrations))) {stop("concentrations must be numeric")}
  if(!all(concentrations >= 0 & is.finite(concentrations))) {
    stop("concentrations must be non-negative real numbers")
  }
  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(length(concentrations) != length(cqs)) {
    stop("concentrations and cqs must be the same length")
  }

  # --- Generate plot

  ggplot2::theme_set(ggplot2::theme_bw())

  plot <- ggplot2::ggplot(esc_data, ggplot2::aes(x = concentrations, y = cqs)) +
    ggplot2::geom_point(na.rm = TRUE, size = 3, alpha = 0.4) +
    ggplot2::scale_x_log10() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank())+
    ggplot2::labs(x = 'Concentration', y = 'Cq')

  return(plot)
}



#' Plot Fitted ESC Model
#'
#' Adds fitted median and probability interval of ESC model to plot of Cq vs
#' Concentration
#'
#' @param model esc object represeting fitted model
#'
#' @return ggplot object representing plot
#' @export
#'
#' @examples
plot_esc_model <- function(model) {

  # --- Input checks
  if(class(model) != "esc") {stop("model is not an esc object")}


  # --- Cosmetics
  col.ci  = 'steelblue2'
  col.med = 'steelblue3'
  alpha.ci = 0.2
  size.med = 1

  # --- Generate plot
  g.data = plot_esc_data(model$data)

  g.model = g.data +
    ggplot2::geom_line(ggplot2::aes(y = model$cq_quantiles$median), color=col.med,
                       linewidth = size.med)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin = model$cq_quantiles$low,
                                      ymax = model$cq_quantiles$high),
                         alpha = alpha.ci, fill = col.ci)+
    ggplot2::labs(title = "ESC model fit",
                  subtitle = paste0('Median and ',
                                    round(model$cq_quantiles_level*100),'%CI'))

  return(g.model)
}
