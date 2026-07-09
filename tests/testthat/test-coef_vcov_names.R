test_that("coef and vcov methods are dimensionally consistent", {
  skip_if_not_installed("VGAM")

  data(pneumo, package = "VGAM")
  pneumo <- transform(pneumo, let = log(exposure.time))

  fit <- tsco(
    # Columns are deliberately not in ordinal order. Supplying `levels`
    # should cause tsco() to reorder grouped counts before fitting.
    cbind(normal, mild, severe) ~ let,
    data = pneumo,
    levels = c("mild", "normal", "severe"),
    cutoff_level = "severe",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  b <- coef(fit)
  V <- vcov(fit)

  expect_equal(names(b), rownames(V))
  expect_equal(names(b), colnames(V))
  expect_equal(nrow(V), length(b))
  expect_equal(ncol(V), length(b))
})
