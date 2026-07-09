test_that("predict.tsco returns valid probabilities", {
  skip_if_not_installed("VGAM")

  data(pneumo, package = "VGAM")
  pneumo <- transform(pneumo, let = log(exposure.time))

  fit <- tsco(
    cbind(normal, mild, severe) ~ let,
    data = pneumo,
    levels = c("mild", "normal", "severe"),
    cutoff_level = "severe",
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
  expect_equal(rowSums(p), rep(1, nrow(p)), tolerance = 1e-8)
})
