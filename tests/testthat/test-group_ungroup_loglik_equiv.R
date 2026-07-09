test_that("grouped and ungrouped TSCO log-likelihoods agree", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")

  data(pneumo, package = "VGAM")
  pneumo <- transform(pneumo, let = log(exposure.time))

  pneumo_ungrouped <- pneumo |>
    tidyr::pivot_longer(
      cols = c(normal, mild, severe),
      names_to = "y",
      values_to = "n"
    ) |>
    tidyr::uncount(weights = n) |>
    dplyr::mutate(
      y = ordered(y, levels = c("mild", "normal", "severe"))
    )

  fit_grouped <- tsco(
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

  fit_ungrouped <- tsco(
    y ~ let,
    data = pneumo_ungrouped,
    levels = c("mild", "normal", "severe"),
    cutoff_level = "severe",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_equal(
    fit_grouped$logLik_stage1,
    fit_ungrouped$logLik_stage1,
    tolerance = 1e-7
  )

  expect_equal(
    fit_grouped$logLik_stage2,
    fit_ungrouped$logLik_stage2,
    tolerance = 1e-7
  )

  expect_equal(
    fit_grouped$logLik,
    fit_ungrouped$logLik,
    tolerance = 1e-7
  )

  expect_equal(
    as.numeric(logLik(fit_grouped)),
    as.numeric(logLik(fit_ungrouped)),
    tolerance = 1e-7
  )

  expect_equal(
    as.numeric(logLik(fit_grouped)),
    fit_grouped$logLik,
    tolerance = 1e-12
  )
})
