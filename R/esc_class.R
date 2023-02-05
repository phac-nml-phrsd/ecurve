new_esc <- function(intercept, slope, sigma) {
  res <- list(intercept = intercept, slope = slope, eff = exp(-1/slope) - 1,
              sigma = sigma)
  structure(res, class = "esc")
}
