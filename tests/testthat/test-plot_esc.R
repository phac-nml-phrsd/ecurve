test_that("data plootting input checks work", {
  df <- data.frame(concentrations = c(1, 2), cps = c(1, 2))
  expect_error(plot_esc_data(df), "must contain cqs column")
  df <- data.frame(concentraions = c(1, 2), cqs = c(1, 2))
  expect_error(plot_esc_data(df), "must contain concentrations column")
  df <- data.frame(concentrations = c(1, 2), cqs = c(1, "a"))
  expect_error(plot_esc_data(df), "cqs must be numeric")
  df <- data.frame(concentrations = c(1, "a"), cqs = c(1, 2))
  expect_error(plot_esc_data(df), "concentrations must be numeric")
  df <- data.frame(concentrations = c(1, 2), cqs = c(1, -1))
  expect_error(plot_esc_data(df), "cqs must be non-negative")
  df <- data.frame(concentrations = c(1, -1), cqs = c(1, 2))
  expect_error(plot_esc_data(df), "concentrations must be non-negative")
})

test_that("model plotting input checks work", {
  concs <- 2^rep(seq(from = -2, to = 10), each = 3)
  cq1 <- 38
  sigma <- 0.5
  eff <- 0.97
  cqs <- sim_cqs(concs, cq1 = cq1, eff = eff, sigma = sigma)
  df <- data.frame(concentrations = concs, cqs = cqs)
  model <- esc_mle(df)
  expect_error(plot_esc_model(5), "model is not an esc object")
  expect_error(plot_esc_model(model, PI = "a"), "PI must be numeric")
  expect_error(plot_esc_model(model, PI = 1.1), "PI must be between 0 and 1")
  expect_error(plot_esc_model(model, approximate = "a"), "approximate must be logical")
})

