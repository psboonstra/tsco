test_that("predict.tsco works with newdata", {
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

  nd <- data.frame(let = seq(min(pneumo$let), max(pneumo$let), length.out = 5))

  p <- predict(fit, newdata = nd, type = "prob")

  expect_equal(nrow(p), nrow(nd))
  expect_equal(colnames(p), fit$levels)
  expect_equal(rowSums(p), rep(1, nrow(nd)), tolerance = 1e-8)
})
