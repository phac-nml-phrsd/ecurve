esc_probdens <- function(concentrations, cqs, intercept, slope, sigma) {
  N0starts <- pmax(qpois(1E-15, concentrations), 1)
  N0ends <- qpois(1E-15, concentrations, lower.tail = FALSE)
  N0min <- min(N0starts)
  N0max <- max(N0ends)
  N0s <- N0min:N0max
  mean_cqs <- intercept + slope * log(N0s)
  probdens <- function(params) {
    sum(dpois(N0s[params[3]:params[4]], params[1]) *
          dnorm(params[2], mean = mean_cqs[params[3]:params[4]], sd = sigma))
  }
  apply(cbind(concentrations, cqs, N0starts - N0min + 1, N0ends - N0min + 1),
        MARGIN = 1, FUN = probdens)
}

esc_log_likelihood <- function(params, concentrations, cqs) {
  if(params[3] < 0) {return(Inf)}
  -sum(log(esc_probdens(concentrations, cqs, params[1], params[2], params[3])))
}

#' Fit ESC Model Using MLE
#'
#' @param concentrations Numeric vector of known sample concentrations
#' @param cqs Numeric vector of corresponding Cq values, with non-detects coded
#' as NaN
#'
#' @return esc object representing fitted model
#' @export
#'
#' @examples
esc_mle <- function(concentrations, cqs) {
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
  detects <- !is.nan(cqs)
  concentrations <- concentrations[detects]
  cqs <- cqs[detects]
  naive_sc <- lm(cqs ~ log(concentrations))
  init <- c(unname(coef(naive_sc)), sigma(naive_sc))
  res <- nlm(esc_log_likelihood, init, concentrations = concentrations, cqs = cqs)
  new_esc(intercept = res$estimate[1], slope = res$estimate[2], sigma = res$estimate[3])
}
