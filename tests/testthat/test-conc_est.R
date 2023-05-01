test_that("concentration MLE works", {
  intercept <- 38
  E <- 0.97
  sigma <- 0.2
  model <- new_esc(intercept = intercept, slope = -1/log10(1 + E), sigma = sigma)
  conc <- 20
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  expect_true(abs(conc_mle(data, model) - conc)/conc < 0.5)
  conc <- 500
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  expect_true(abs(conc_mle(data, model) - conc)/conc < 0.5)
  conc <- 1.5
  data <- sim_cqs(rep(conc, 10), eff = E, cq1 = intercept,
                  sigma = sigma)
  expect_true(abs(conc_mle(data, model) - conc)/conc < 1)
})

test_that("concentration interval estimation works", {
  intercept <- 38
  E <- 0.97
  sigma <- 0.2
  model <- new_esc(intercept = intercept, slope = -1/log10(1 + E), sigma = sigma)
  conc <- 20
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_interval(data, model)
  expect_s3_class(int, "conc_int")
  expect_true(abs(int$interval$mle - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)
  int1 <- conc_interval(data, model, level = 0.97)
  expect_equal(int1$interval$mle, int$interval$mle)
  expect_true(int1$interval$upper > int$interval$upper)
  expect_true(int1$interval$lower < int$interval$lower)
  conc <- 500
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_interval(data, model)
  expect_s3_class(int, "conc_int")
  expect_true(abs(int$interval$mle - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)
  conc <- 1.5
  data <- sim_cqs(rep(conc, 10), eff = E, cq1 = intercept,
                  sigma = sigma)
  int <- conc_interval(data, model)
  expect_s3_class(int, "conc_int")
  expect_true(abs(int$interval$mle - conc)/conc < 1)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)
})

test_that("concentration MLE works when approximation is turned off", {
  intercept <- 38
  E <- 0.97
  sigma <- 0.2
  model <- new_esc(intercept = intercept, slope = -1/log10(1 + E), sigma = sigma)
  conc <- 20
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  expect_true(abs(conc_mle(data, model, approximate = FALSE) - conc)/conc < 0.5)
  conc <- 500
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  expect_true(abs(conc_mle(data, model, approximate = FALSE) - conc)/conc < 0.5)
  conc <- 1.5
  data <- sim_cqs(rep(conc, 10), eff = E, cq1 = intercept,
                  sigma = sigma)
  expect_true(abs(conc_mle(data, model, approximate = FALSE) - conc)/conc < 1)
})

test_that("concentration interval estimation works when approximation is turned off", {
  intercept <- 38
  E <- 0.97
  sigma <- 0.2
  model <- new_esc(intercept = intercept, slope = -1/log10(1 + E), sigma = sigma)
  conc <- 20
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_interval(data, model, approximate = FALSE)
  expect_s3_class(int, "conc_int")
  expect_true(abs(int$interval$mle - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)
  int1 <- conc_interval(data, model, level = 0.97, approximate = FALSE)
  expect_equal(int1$interval$mle, int$interval$mle)
  expect_true(int1$interval$upper > int$interval$upper)
  expect_true(int1$interval$lower < int$interval$lower)
  conc <- 500
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_interval(data, model, approximate = FALSE)
  expect_s3_class(int, "conc_int")
  expect_true(abs(int$interval$mle - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)

  conc <- 1.5
  data <- sim_cqs(rep(conc, 10), eff = E, cq1 = intercept,
                  sigma = sigma)
  int <- conc_interval(data, model, approximate = FALSE)
  expect_s3_class(int, "conc_int")
  expect_true(abs(int$interval$mle - conc)/conc < 1)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)
})

test_that("Concentration MLE works if all inputs are NaN",{
  model <- new_esc(intercept = 38, slope = -1/log10(1.95), sigma = 0.5)
  expect_equal(conc_mle(c(NaN, NaN, NaN), model), 0)
})

test_that("Concentration Interval Estimation works if all inputs are NaN",{
  model <- new_esc(intercept = 38, slope = -1/log10(1.95), sigma = 0.5)
  int <- conc_interval(c(NaN, NaN, NaN), model)
  expect_s3_class(int, "conc_int")
  expect_equal(int$interval$lower, 0)
  expect_equal(int$interval$mle, 0)
  expect_equal(1 - exp(-3 * int$interval$upper), 0.95)
})

test_that("conc_mle input checks work", {
  model <- new_esc(intercept = 38, slope = -1/log10(1.95), sigma = 0.5)
  expect_error(conc_mle(c(1, "a", 1), model), "cqs must be numeric")
  expect_error(conc_mle(c(1, -1, 1), model), "cqs must be non-negative")
  expect_error(conc_mle(c(1, 1, 1), 5), "model is not an esc object")
  expect_error(conc_mle(c(1, 1, 1), model, approximate = "a"),
               "approximate must be logical")
})

test_that("conc_interval input checks work", {
  model <- new_esc(intercept = 38, slope = -1/log10(1.95), sigma = 0.5)
  expect_error(conc_interval(c(1, "a", 1), model), "cqs must be numeric")
  expect_error(conc_interval(c(1, -1, 1), model), "cqs must be non-negative")
  expect_error(conc_interval(c(1, 1, 1), 5), "model is not an esc object")
  expect_error(conc_interval(c(1, 1), model, level = "a"), "level must be numeric")
  expect_error(conc_interval(c(1, 1), model, level = 1.2),
               "level must be between 0 and 1")
  expect_error(conc_interval(c(1, 1), model, approximate = "a"),
               "approximate must be logical")
})

test_that("multi_interval works", {
  intercept <- 38
  E <- 0.97
  sigma <- 0.2
  model <- new_esc(intercept, -1/log10(1 + E), sigma)
  conc <- 20
  cqs <- sim_cqs(rep(conc, 5), E, intercept, sigma)
  df <- data.frame(sample = c(rep(1, 5), rep(2, 3)), cqs = c(cqs, rep(NaN, 3)))
  res <- multi_interval(df, model)
  expect_equal(dim(res), c(2, 5))
  expect_setequal(res$sample, c(1, 2))
  ind1 <- which(res$sample == 1)
  ind2 <- which(res$sample == 2)
  expect_true(abs(res$mle[ind1] - conc)/conc < 0.5)
  expect_true(res$lower[ind1] < res$mle[ind1])
  expect_true(res$upper[ind1] > res$mle[ind1])
  expect_equal(res$mle[ind2], 0)
  expect_equal(res$lower[ind2], 0)
  expect_equal(1 - exp(-3 * res$upper[ind2]), 0.95)
  df <- data.frame(sample = c(rep("a", 5), rep("b", 3)), cqs = c(cqs, rep(NaN, 3)))
  res <- multi_interval(df, model, level = 0.90, approximate = FALSE)
  expect_equal(dim(res), c(2, 5))
  expect_setequal(res$sample, c("a", "b"))
  ind1 <- which(res$sample == "a")
  ind2 <- which(res$sample == "b")
  expect_true(abs(res$mle[ind1] - conc)/conc < 0.5)
  expect_true(res$lower[ind1] < res$mle[ind1])
  expect_true(res$upper[ind1] > res$mle[ind1])
  expect_equal(res$mle[ind2], 0)
  expect_equal(res$lower[ind2], 0)
  expect_equal(1 - exp(-3 * res$upper[ind2]), 0.90)
})

test_that("multi_interval input checks work", {
  model <- new_esc(38, -1/log(1.97), 0.5)
  df <- data.frame(index = c(1, 1, 2, 2), cqs = c(1, 2, 1, 2))
  expect_error(multi_interval(df, model), "cq_data must contain sample column")
  df <- data.frame(sample = c(1, 1, 2, 2), cps = c(1, 2, 1, 2))
  expect_error(multi_interval(df, model), "cq_data must contain cqs column")
  df <- data.frame(sample = c(1, 1, 2, 2), cqs = c(1, 2, "a", 2))
  expect_error(multi_interval(df, model), "cqs must be numeric")
  df <- data.frame(sample = c(1, 1, 2, 2), cqs = c(1, 2, -5, 2))
  expect_error(multi_interval(df, model), "cqs must be non-negative")
  df <- data.frame(sample = c(1, 1, 2, 2), cqs = c(1, 2, 1, 2))
  expect_error(multi_interval(df, 5), "model is not an esc object")
  expect_error(multi_interval(df, model, level = "a"), "level must be numeric")
  expect_error(multi_interval(df, model, level = 1.5),
               "level must be between 0 and 1")
  expect_error(multi_interval(df, model, approximate = "a"),
               "approximate must be logical")
})

test_that("concentration mcmc estimation works", {
  skip_if_not_installed("runjags")

  intercept <- 38
  E         <- 0.97
  sigma     <- 0.2
  model <- new_esc(intercept = intercept, slope = -1/log10(1 + E), sigma = sigma)

  conc <- 20
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_mcmc(data, model)
  expect_s3_class(int$mcmc_samples, "runjags")
  expect_true(abs(int$interval$median - conc)/conc < 0.5)
  expect_true(abs(int$interval$mean - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$median)
  expect_true(int$interval$lower < int$interval$median)
  int1 <- conc_mcmc(data, model, level = 0.99)
  expect_true(abs(int1$interval$median - conc)/conc < 0.5)
  expect_true(abs(int1$interval$mean - conc)/conc < 0.5)
  expect_true(int1$interval$upper > int$interval$upper)
  expect_true(int1$interval$lower < int$interval$lower)

  conc <- 500
  data <- sim_cqs(rep(conc, 5), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_mcmc(data, model)
  expect_s3_class(int$mcmc_samples, "runjags")
  expect_true(abs(int$interval$median - conc)/conc < 0.5)
  expect_true(abs(int$interval$mean - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$median)
  expect_true(int$interval$lower < int$interval$median)

  conc <- 1.5
  data <- sim_cqs(rep(conc, 10), eff = E, cq1 = intercept,
                  sigma = sigma)
  int <- conc_mcmc(data, model)
  expect_s3_class(int$mcmc_samples, "runjags")
  expect_true(abs(int$interval$median - conc)/conc < 1)
  expect_true(abs(int$interval$mean - conc)/conc < 1)
  expect_true(int$interval$upper > int$interval$median)
  expect_true(int$interval$lower < int$interval$median)
})

test_that("conc_mcmc input checks work", {
  skip_if_not_installed("runjags")
  model <- new_esc(intercept = 38, slope = -1/log10(1.95), sigma = 0.5)
  expect_error(conc_mcmc(c(1, "a", 1), model), "cqs must be numeric")
  expect_error(conc_mcmc(c(1, -1, 1), model), "cqs must be non-negative")
  expect_error(conc_mcmc(c(1, 1, 1), 5), "model is not an esc object")
  expect_error(conc_mcmc(c(1, 1), model, level = "a"), "level must be numeric")
  expect_error(conc_mcmc(c(1, 1), model, level = 1.2),
               "level must be between 0 and 1")
})

test_that("multi_conc_mcmc works", {
  skip_if_not_installed("runjags")
  intercept <- 38
  E <- 0.97
  sigma <- 0.2
  model <- new_esc(intercept, -1/log10(1 + E), sigma)
  conc <- 20
  cqs <- sim_cqs(rep(conc, 5), E, intercept, sigma)
  df <- data.frame(sample = c(rep(1, 5), rep(2, 3)), cqs = c(cqs, rep(NaN, 3)))
  res <- multi_conc_mcmc(df, model)
  expect_s3_class(res$mcmc_samples, "runjags")
  res <- res$intervals
  expect_equal(dim(res), c(2, 5))
  expect_setequal(res$sample, c(1, 2))
  ind1 <- which(res$sample == 1)
  ind2 <- which(res$sample == 2)
  expect_true(abs(res$median[ind1] - conc)/conc < 0.5)
  expect_true(abs(res$mean[ind1] - conc)/conc < 0.5)
  expect_true(res$lower[ind1] < res$median[ind1])
  expect_true(res$upper[ind1] > res$median[ind1])
  expect_true(res$lower[ind2] < 1e-4)
  expect_true(res$median[ind2] < 1)
  expect_true(res$mean[ind2] < 1)
  expect_true(res$upper[ind2] < 2)
  df <- data.frame(sample = c(rep("a", 5), rep("b", 3)), cqs = c(cqs, rep(NaN, 3)))
  res <- multi_conc_mcmc(df, model, level = 0.90)
  expect_s3_class(res$mcmc_samples, "runjags")
  res <- res$intervals
  expect_equal(dim(res), c(2, 5))
  expect_setequal(res$sample, c("a", "b"))
  ind1 <- which(res$sample == "a")
  ind2 <- which(res$sample == "b")
  expect_true(abs(res$median[ind1] - conc)/conc < 0.5)
  expect_true(abs(res$mean[ind1] - conc)/conc < 0.5)
  expect_true(res$lower[ind1] < res$median[ind1])
  expect_true(res$upper[ind1] > res$median[ind1])
  expect_true(res$lower[ind2] < 1e-4)
  expect_true(res$median[ind2] < 1)
  expect_true(res$mean[ind2] < 1)
  expect_true(res$upper[ind2] < 2)
})

test_that("multi_conc_mcmc input checks work", {
  skip_if_not_installed("runjags")
  model <- new_esc(38, -1/log(1.97), 0.5)
  df <- data.frame(index = c(1, 1, 2, 2), cqs = c(1, 2, 1, 2))
  expect_error(multi_conc_mcmc(df, model), "cq_data must contain sample column")
  df <- data.frame(sample = c(1, 1, 2, 2), cps = c(1, 2, 1, 2))
  expect_error(multi_conc_mcmc(df, model), "cq_data must contain cqs column")
  df <- data.frame(sample = c(1, 1, 2, 2), cqs = c(1, 2, "a", 2))
  expect_error(multi_conc_mcmc(df, model), "cqs must be numeric")
  df <- data.frame(sample = c(1, 1, 2, 2), cqs = c(1, 2, -5, 2))
  expect_error(multi_conc_mcmc(df, model), "cqs must be non-negative")
  df <- data.frame(sample = c(1, 1, 2, 2), cqs = c(1, 2, 1, 2))
  expect_error(multi_conc_mcmc(df, 5), "model is not an esc object")
  expect_error(multi_conc_mcmc(df, model, level = "a"), "level must be numeric")
  expect_error(multi_conc_mcmc(df, model, level = 1.5),
               "level must be between 0 and 1")
})
