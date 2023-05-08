#' @title Plot Standard Curve Calibration Data
#'
#' @description Given a data frame containing concentrations
#' and corresponding Cq values, plots Cq values against
#' log-transformed concentrations.
#'
#' @param esc_data Data frame containing data to bt used to fit the model. Must
#' contain a column named \code{concentrations} with known sample concentrations,
#'  and a column named \code{cqs} with corresponding Cq values.
#'  Non-detects should be encoded by a Cq value of \code{NaN}.
#' @param xlimits Optional numeric vector of length 2. If provided, passed to
#' \code{ggplot2} to specify x-axis limits of plot
#'
#' @return A \code{ggplot} object representing plot.
#'
#' @export
#'
#' @examples
#' esc_data = data.frame(
#'     concentrations = c(1,1,10, 10, 100, 500, 500),
#'     cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
#' plot_esc_data(esc_data)
#'
plot_esc_data <- function(esc_data, xlimits = NULL) {

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

  # --- Generate plot

  ggplot2::theme_set(ggplot2::theme_bw())

  plot <- ggplot2::ggplot(esc_data, ggplot2::aes(x = concentrations, y = cqs)) +
    ggplot2::geom_point(na.rm = TRUE, size = 3, alpha = 0.4) +
    ggplot2::scale_x_continuous(trans = "log10", limits = xlimits) +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank())+
    ggplot2::labs(x = 'Concentration (gc/rxn)', y = 'Cq')

  return(plot)
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
cq_quantile <- function(alpha,
                        concentrations,
                        intercept,
                        slope,
                        sigma) {
  N0mins <- pmax(qpois(1E-15, concentrations), 1)
  N0maxes <- pmax(qpois(1E-15, concentrations, lower.tail = FALSE), N0mins)
  N0s <- unique(unlist(mapply(FUN = seq, sort(N0mins), sort(N0maxes))))
  means <- intercept + slope * log(N0s)
  N0starts <- match(N0mins, N0s)
  N0ends <- match(N0maxes, N0s)

  quant <- function(conc, N0start, N0end) {
    stats::uniroot(function(cq) {
      sum(pnorm(cq, mean = means[N0start:N0end], sd = sigma) *
            dpois(N0s[N0start:N0end], conc)) - alpha * (1 - exp(-conc))},
      lower = qnorm(alpha, mean = means[N0end], sd = sigma),
      upper = qnorm(alpha, mean = means[N0start], sd = sigma))$root
  }
  res = mapply(quant, conc = concentrations, N0start = N0starts, N0end = N0ends)
  return(res)
}

#' Faster Estimation of Theoretical Quantiles for Cq Values
#'
#' Given ESC model parameters, quantile level alpha, and concentration conc,
#' estimates that alpha-th quantile of the theoretical distribution of Cq values
#' at the specified concentration. For efficiency, ensures that a maximum of 100
#' N0 values are considered when the computation is performed. N0 values used for
#' computation are equally spaced out between the computed lower and upper bounds,
#' and are used to estimate contribution of nearby N0 values not included in the
#' computation.
#'
#' @param conc numeric concentration at which to evaluate specified
#' quantile
#' @param alpha quantile level
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return estimate of specified quantile
cq_quantile_est <- function(alpha, conc, intercept, slope, sigma) {
  #determine N0 values at which to perform evaluation
  N0start <- max(qpois(1e-15, conc), 1)
  N0end <- max(qpois(1e-15, conc, lower.tail = FALSE), N0start + 1)
  granularity <- ceiling((N0end - N0start)/100)
  N0s <- seq(N0start, N0end, by = granularity)

  #compute quantile
  stats::uniroot(function(cq) {
    sum(pnorm(cq, mean = intercept + slope * log(N0s), sd = sigma) *
          dpois(N0s, conc) * granularity) - alpha * (1 - exp(-conc))},
    upper = qnorm(alpha, mean = intercept + slope * log(N0start), sd = sigma),
    lower = qnorm(alpha, mean = intercept + slope * log(N0end), sd = sigma))$root
}


#' @title Plot Fitted ESC Model
#'
#' @description Plots the fitted median and probability interval of
#' an ESC model Cq vs concentration.
#'
#' @param model \code{esc} object representing fitted model.
#' @param PI Numeric. Width of the probability interval (must be between 0 and 1).
#' @param title String. The title of the plot. If unspecified, defaults to "ESC
#' model fit"
#' @param xlimits Optional numeric vector of length 2. If provided, passed to
#' \code{ggplot2} to specify x-axis limits of plot
#' @param approximate Logical. If \code{TRUE} (the default), a faster but potentially
#' less accurate approximation for the likelihood function will be used at high
#' concentrations.
#'
#' @return A \code{ggplot} object representing plot.
#' @export
#'
#' @examples
#' esc_data = data.frame(
#'     concentrations = c(1,1,10, 10, 100, 500, 500),
#'     cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
#' mod = esc_mle(esc_data)
#' plot_esc_model(mod)
#'
#'
plot_esc_model <- function(model, PI = 0.95, title = "ESC model fit",
                           xlimits = NULL, approximate = TRUE) {

  # --- Input checks
  if(!inherits(model, 'esc')) {stop("model is not an esc object")}

  if(is.null(model$data)){
    message('The ESC model object input does not contain `data`. The function `plot_esc_model()` can only plot ESC object with data. Returning no plot.')
   return(NULL)
  }

  if(!is.numeric(PI)) {stop("PI must be numeric")}
  if(PI > 1 | PI < 0) {stop("PI must be between 0 and 1")}

  if(!is.logical(approximate)) {stop("approximate must be logical")}

  if(!is.character(title) | !length(title) == 1) {stop("title must be a string")}

  # --- Cosmetics

  col.ci  = 'steelblue2'
  col.med = 'steelblue3'
  alpha.ci = 0.2
  size.med = 1

  # --- Quantiles for Cq values

  alphas <- c(0.5 - PI/2, 0.5, 0.5 + PI/2)

  llim <- log(min(model$data$concentrations))
  ulim <- log(max(model$data$concentrations))
  xpsn_factor <- (ulim - llim) * 0.05
  concentrations <- exp(seq(llim - xpsn_factor, ulim + xpsn_factor,
                            length.out = 100))

  cq_quants <- lapply(X = alphas,
                     FUN = ifelse(approximate,
                                  Vectorize(cq_quantile_est, "conc"),
                                  cq_quantile),
                     concentrations,
                     intercept = model$intercept,
                     slope = model$slope / log(10),
                     sigma = model$sigma)

  names(cq_quants) <- c('low', 'median', 'high')
  cq_quants <- as.data.frame(cq_quants)

  # --- Generate plot
  g.data = plot_esc_data(model$data, xlimits)

  g.model = g.data +
    ggplot2::geom_line(data = cq_quants, ggplot2::aes(x = concentrations,
                                                      y = median),
                       color=col.med, linewidth = size.med)+
    ggplot2::geom_ribbon(data = cq_quants, ggplot2::aes(x = concentrations,
                                                        ymin = low, ymax = high,
                                                        y = median),
                         alpha = alpha.ci, fill = col.ci)+
    ggplot2::labs(title = title,
                  subtitle = paste0('Median and ',
                                    round(PI*100),'%PI'))

  return(g.model)
}
