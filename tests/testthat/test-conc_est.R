test_that("concentration MLE works", {
  intercept <- 38
  E <- 0.97
  sigma <- 0.5
  model <- new_esc(intercept = intercept, slope = -1/log(1 + E), sigma = sigma,
                   cq_quantiles = NULL, cq_quantiles_CI = NULL)
  conc <- 20
  data <- sim_cqs(c(conc, conc, conc), eff = E, cq1 = intercept, sigma = sigma)
  expect_true(abs(conc_mle(data, model) - conc)/conc < 0.5)
  conc <- 500
  data <- sim_cqs(c(conc, conc, conc), eff = E, cq1 = intercept, sigma = sigma)
  expect_true(abs(conc_mle(data, model) - conc)/conc < 0.5)
  conc <- 1.5
  data <- sim_cqs(c(conc, conc, conc, conc, conc, conc), eff = E, cq1 = intercept,
                  sigma = sigma)
  expect_true(abs(conc_mle(data, model) - conc)/conc < 0.5)
})

test_that("concentration interval estimation works", {
  intercept <- 38
  E <- 0.97
  sigma <- 0.5
  model <- new_esc(intercept = intercept, slope = -1/log10(1 + E), sigma = sigma)
  conc <- 20
  data <- sim_cqs(c(conc, conc, conc), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_interval(data, model)
  expect_s3_class(int, "conc_int")
  expect_true((int$interval$mle - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)
  int1 <- conc_interval(data, model, level = 0.97)
  expect_equal(int1$interval$mle, int$interval$mle)
  expect_true(int1$interval$upper > int$interval$upper)
  expect_true(int1$interval$lower < int$interval$lower)
  conc <- 500
  data <- sim_cqs(c(conc, conc, conc), eff = E, cq1 = intercept, sigma = sigma)
  int <- conc_interval(data, model)
  expect_s3_class(int, "conc_int")
  expect_true((int$interval$mle - conc)/conc < 0.5)
  expect_true(int$interval$upper > int$interval$mle)
  expect_true(int$interval$lower < int$interval$mle)
  conc <- 1.5
  data <- sim_cqs(c(conc, conc, conc, conc, conc, conc), eff = E, cq1 = intercept,
                  sigma = sigma)
  int <- conc_interval(data, model)
  expect_s3_class(int, "conc_int")
  expect_true((int$interval$mle - conc)/conc < 0.5)
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
