test_that("stage-specific coef and vcov methods agree with underlying VGAM fits", {
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

  expect_equal(
    coef(fit, stage = "stage1"),
    stats::coef(fit$fit_stage1)
  )

  expect_equal(
    coef(fit, stage = "stage2"),
    stats::coef(fit$fit_stage2)
  )

  expect_equal(
    vcov(fit, stage = "stage1"),
    stats::vcov(fit$fit_stage1)
  )

  expect_equal(
    vcov(fit, stage = "stage2"),
    stats::vcov(fit$fit_stage2)
  )
})
