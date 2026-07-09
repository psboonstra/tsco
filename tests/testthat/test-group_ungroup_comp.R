test_that("grouped and ungrouped fits agree when levels reorder grouped columns", {
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
    # outcomes are deliberately not cbinded correctly, to test
    # robustness to column order (as long as levels is provided)
    cbind(normal, mild, severe) ~ let,
    data = pneumo,
    levels = c("mild", "normal", "severe"),
    cutoff_level = "severe",
    stage1 = "po",
    stage2 = "po"
  )

  fit_ungrouped <- tsco(
    y ~ let,
    data = pneumo_ungrouped,
    levels = c("mild", "normal", "severe"),
    cutoff_level = "severe",
    stage1 = "po",
    stage2 = "po"
  )

  expect_equal(
    unname(coef(fit_grouped$fit_stage1)),
    unname(coef(fit_ungrouped$fit_stage1)),
    tolerance = 1e-7
  )

  expect_equal(
    unname(coef(fit_grouped$fit_stage2)),
    unname(coef(fit_ungrouped$fit_stage2)),
    tolerance = 1e-7
  )
})
