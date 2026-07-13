test_that("predict.tsco returns valid wine probabilities", {
  dat <- make_wine_test_data()

  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  p <- predict(fit, type = "prob")

  expect_true(is.matrix(p))
  expect_equal(colnames(p), fit$levels)
  expect_equal(nrow(p), fit$n_groups)
  expect_true(all(p >= -1e-10))
  expect_true(all(p <= 1 + 1e-10))
  expect_equal(rowSums(p), rep(1, nrow(p)), tolerance = 1e-8, ignore_attr = TRUE)
})
