test_that("coef and vcov methods are dimensionally consistent", {
  dat <- make_wine_test_data()

  # This is the manuscript's headline PO|C|MR variant: conditional PO
  # below the cutoff, marginal (saturated) multinomial at/above it.
  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  b <- coef(fit)
  V <- vcov(fit)

  expect_equal(names(b), rownames(V))
  expect_equal(names(b), colnames(V))
  expect_equal(nrow(V), length(b))
  expect_equal(ncol(V), length(b))
})
