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
  plot <- ggplot2::ggplot(data.frame(concentration = concentrations, cq = cqs),
                          ggplot2::aes(x = concentration, y = cq))
  plot <- plot + ggplot2::geom_point(na.rm = TRUE) + ggplot2::scale_x_log10()
  print(plot)
  plot
}

cq_quantile_factory <- function(alpha, intercept, slope, sigma) {
  function(concentrations) {
    N0starts <- pmax(qpois(1E-15, concentrations), 1)
    N0ends <- pmax(qpois(1E-15, concentrations, lower.tail = FALSE), N0starts + 1)
    N0min <- min(N0starts)
    N0max <- max(N0ends)
    N0s <- N0min:N0max
    means <- intercept + slope * log(N0s)
    cq_quantile <- function(params) {
      uniroot(function(cq) {
        sum(pnorm(cq, mean = means[params[2]:params[3]], sd = sigma) *
              dpois(N0s[params[2]:params[3]], params[1])) -
          alpha * (1 - exp(-params[1]))},
                lower = qnorm(alpha, mean = means[params[3]], sd = sigma),
                upper = qnorm(alpha, mean = means[params[2]], sd = sigma))$root
    }
    apply(cbind(concentrations, N0starts - N0min + 1, N0ends - N0min + 1),
          MARGIN = 1, FUN = cq_quantile)
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
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}
  plot <- plot +
    ggplot2::geom_function(fun = cq_quantile_factory(0.5, model$intercept,
                                                     model$slope, model$sigma),
                           color = "red") +
    ggplot2::geom_function(fun = cq_quantile_factory((1 - level)/2, model$intercept,
                                                     model$slope, model$sigma),
                           color = "blue") +
    ggplot2::geom_function(fun = cq_quantile_factory(1 - (1 - level)/2, model$intercept,
                                                     model$slope, model$sigma),
                           color = "blue")
  print(plot)
  plot
}
