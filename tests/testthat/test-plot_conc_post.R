test_that("interval plotting input checks work", {

  esc_data = data.frame(
    concentrations = c(1,1,10, 10, 100, 500, 500),
    cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
  mod = esc_mle(esc_data)
  new.cqs = c( 35, NaN, 36)
  x = conc_interval(cqs = new.cqs, model = mod)
  g1 = plot_conc_post(interval = x, type = 'pdf')
  g2 = plot_conc_post(interval = x, type = 'cdf')


  expect_s3_class(g1, 'ggplot')
  expect_s3_class(g2, 'ggplot')

  expect_error(plot_conc_post(5, "pdf"), "interval is not a conc_int object")
  expect_error(plot_conc_post(interval = x, "vdf"), "invalid type")
  expect_error(plot_conc_post(interval = x, type = 'cdf', title = matrix(1:9)))
})
