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

  #Input checks
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

  #Generate plot
  plot <- ggplot2::ggplot(data.frame(concentration = concentrations, cq = cqs),
                          ggplot2::aes(x = concentration, y = cq))
  plot <- plot + ggplot2::geom_point(na.rm = TRUE) + ggplot2::scale_x_log10()

  plot
}

#' Calculate Theoretical Quantiles for Cq Values
#'
#' Given ESC model parameters, quantile level alpha, and vector of concentrations
#' evaluates that alpha-th quantile of the theoretical distribution of cq values
#' at each of the specified concentrations
#'
#' @param concentrations numeric vector of concentrations at which to evaluate
#' specified quantile
#' @param alpha quantile level
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return Numeric vector containing the calculated quantiles
cq_quantile <- function(concentrations, alpha, intercept, slope, sigma) {
  N0starts <- pmax(qpois(1E-15, concentrations), 1)
  N0ends <- pmax(qpois(1E-15, concentrations, lower.tail = FALSE), N0starts + 1)
  N0min <- min(N0starts)
  N0max <- max(N0ends)
  N0s <- N0min:N0max
  means <- intercept + slope * log(N0s)
  quant <- function(conc, N0start, N0end) {
    uniroot(function(cq) {
      sum(pnorm(cq, mean = means[N0start:N0end], sd = sigma) *
            dpois(N0s[N0start:N0end], conc)) - alpha * (1 - exp(-conc))},
              lower = qnorm(alpha, mean = means[N0end], sd = sigma),
              upper = qnorm(alpha, mean = means[N0start], sd = sigma))$root
  }
  mapply(quant, conc = concentrations, N0start = N0starts - N0min + 1,
          N0end = N0ends - N0min + 1)
}

#' Plot Fitted ESC Model
#'
#' Adds fitted median and probability interval of ESC model to plot of Cq vs
#' Concentration
#'
#' @param model esc object represeting fitted model
#' @param plot ggplot object to which fitted model will be added. Optional,
#' if none is provided, a plot of the data used to fit the model will be
#' automatically generated using [plot_esc_data()].
#' @param level level of probability interval, defaults to 0.95
#'
#' @return ggplot object representing plot
#' @export
#'
#' @examples
plot_esc_model <- function(model, plot, level = 0.95) {

  #Input checks
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}

  if(missing(plot)) {plot <- plot_esc_data(model$data)}

  #Modify plot
  plot <- plot +
    ggplot2::geom_function(fun = cq_quantile,
                           args = list(alpha = 0.5,
                                       intercept = model$intercept,
                                       slope = model$slope,
                                       sigma = model$sigma),
                           ggplot2::aes(color = "Median")) +
    ggplot2::geom_function(fun = cq_quantile,
                           args = list(alpha = (1 - level)/2,
                                       intercept = model$intercept,
                                       slope = model$slope,
                                       sigma = model$sigma),
                           ggplot2::aes(color = "Interval")) +
    ggplot2::geom_function(fun = cq_quantile,
                           args = list(alpha = 1 - (1 - level)/2,
                                       intercept = model$intercept,
                                       slope = model$slope,
                                       sigma = model$sigma),
                           ggplot2::aes(color = "Interval")) +
    ggplot2::scale_color_manual(name = "Legend",
                                values = c("Median" = "red", "Interval" = "blue"))

  plot
}
