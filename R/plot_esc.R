#' Plot Standard Curve Calibration Data
#'
#' Given vectors of concentrations and corresponding Cq values, plots Cq values
#' agains log-transformed concentrations
#'
#' @param concentrations Numeric vector of known sample concentrations
#' @param cqs Numeric vector of corresponding Cq values, with non-detects coded
#' as NaN
#'
#' @return ggplot object representing plot
#' @export
#'
#' @examples
plot_esc_data <- function(concentrations, cqs) {

  #Input checks
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

  #Display and return plot
  print(plot)
  plot
}

#' Factory Function for Generating Cq Quantile Function
#'
#' Given ESC model parameters and quantile level alpha, returns function giving
#' corresponding theoretical quantile of Cq values as a function of concentration
#'
#' @param alpha quantile level
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return Function taking single input concentration and returning corresponding
#' Cq qunatile
cq_quantile_factory <- function(alpha, intercept, slope, sigma) {
  function(concentrations) {
    N0starts <- pmax(qpois(1E-15, concentrations), 1)
    N0ends <- pmax(qpois(1E-15, concentrations, lower.tail = FALSE), N0starts + 1)
    N0min <- min(N0starts)
    N0max <- max(N0ends)
    N0s <- N0min:N0max
    means <- intercept + slope * log(N0s)
    cq_quantile <- function(conc, N0start, N0end) {
      uniroot(function(cq) {
        sum(pnorm(cq, mean = means[N0start:N0end], sd = sigma) *
              dpois(N0s[N0start:N0end], conc)) - alpha * (1 - exp(-conc))},
                lower = qnorm(alpha, mean = means[N0end], sd = sigma),
                upper = qnorm(alpha, mean = means[N0start], sd = sigma))$root
    }
    mapply(cq_quantile, conc = concentrations, N0start = N0starts - N0min + 1,
           N0end = N0ends - N0min + 1)
  }
}

#' Plot Fitted ESC Model
#'
#' Adds fitted median and probability interval of ESC model to plot of Cq vs
#' Concentration
#'
#' @param plot ggplot object to which fitted model will be added
#' @param model esc object represeting fitted model
#' @param level level of probability interval, defaults to 0.95
#'
#' @return ggplot object representing plot
#' @export
#'
#' @examples
plot_esc_model <- function(plot, model, level = 0.95) {

  #Input checks
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}

  #Modify plot
  plot <- plot +
    ggplot2::geom_function(fun = cq_quantile_factory(0.5, model$intercept,
                                                     model$slope, model$sigma),
                           ggplot2::aes(color = "Median")) +
    ggplot2::geom_function(fun = cq_quantile_factory((1 - level)/2, model$intercept,
                                                     model$slope, model$sigma),
                           ggplot2::aes(color = "Interval")) +
    ggplot2::geom_function(fun = cq_quantile_factory(1 - (1 - level)/2, model$intercept,
                                                     model$slope, model$sigma),
                           ggplot2::aes(color = "Interval")) +
    ggplot2::scale_color_manual(name = "Legend",
                                values = c("Median" = "red", "Interval" = "blue"))

  #Display and return plot
  print(plot)
  plot
}
