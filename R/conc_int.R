new_conc_int <- function(mle, interval, grid, pdf, cdf) {
  norm_factor <- cdf[length(cdf)]
  structure(list(
    interval = list(lower = interval[1], mle = mle, upper = interval[2]),
    distribution = data.frame(concentration = grid, pdf = pdf/norm_factor,
                              cdf = cdf/norm_factor)
  ), class = "conc_int")
}
