test_that("grouped and individual wine fits agree", {
  dat <- make_wine_test_data()

  fit_grouped <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  fit_ind <- tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_equal(
    unname(coef(fit_grouped)),
    unname(coef(fit_ind)),
    tolerance = 1e-5
  )

  expect_equal(
    as.numeric(logLik(fit_grouped)),
    as.numeric(logLik(fit_ind)),
    tolerance = 1e-7
  )

  expect_equal(
    as.numeric(AIC(fit_grouped)),
    as.numeric(AIC(fit_ind)),
    tolerance = 1e-7
  )

  expect_equal(
    as.numeric(BIC(fit_grouped)),
    as.numeric(BIC(fit_ind)),
    tolerance = 1e-7
  )
})
