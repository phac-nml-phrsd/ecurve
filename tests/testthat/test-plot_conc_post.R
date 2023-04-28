test_that("interval plotting input checks work", {
  model <- new_esc(intercept = 38, slope = -1/log10(1.95), sigma = 0.5)
  int <- conc_interval(c(34, 34.1, 34.2), model)
  expect_error(plot_conc_post(5, "pdf"), "interval is not a conc_int object")
  expect_error(plot_conc_post(int, "vdf"), "invalid type")
})
