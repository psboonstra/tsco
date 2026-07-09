test_that("AIC and BIC are consistent across methods and data representations", {
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

  s_grouped <- summary(fit_grouped)
  s_ungrouped <- summary(fit_ungrouped)

  expect_equal(
    as.numeric(stats::AIC(fit_grouped)),
    s_grouped$total_fit$AIC,
    tolerance = 1e-8
  )

  expect_equal(
    as.numeric(stats::BIC(fit_grouped)),
    s_grouped$total_fit$BIC,
    tolerance = 1e-8
  )

  expect_equal(
    as.numeric(stats::AIC(fit_grouped)),
    as.numeric(stats::AIC(fit_ungrouped)),
    tolerance = 1e-7
  )

  expect_equal(
    as.numeric(stats::BIC(fit_grouped)),
    as.numeric(stats::BIC(fit_ungrouped)),
    tolerance = 1e-7
  )

  expect_equal(
    s_grouped$total_fit$AIC,
    s_ungrouped$total_fit$AIC,
    tolerance = 1e-7
  )

  expect_equal(
    s_grouped$total_fit$BIC,
    s_ungrouped$total_fit$BIC,
    tolerance = 1e-7
  )
})
