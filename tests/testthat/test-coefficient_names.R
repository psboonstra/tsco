test_that("multinomial coefficient names use outcome labels", {
  dat <- make_wine_test_data()

  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "multinomial",
    warn_degenerate = FALSE
  )

  b2 <- coef(fit, stage = "stage2")
  V2 <- vcov(fit, stage = "stage2")
  s <- summary(fit, joint_test = "none")

  raw_names <- names(stats::coef(fit$fit_stage2))
  expected_categories <- fit$collapsed_levels[-1L]
  expected_names <- raw_names

  for (j in seq_along(expected_categories)) {
    expected_names <- sub(
      paste0(":", j, "$"),
      paste0(":", expected_categories[j]),
      expected_names
    )
  }

  expect_identical(names(b2), expected_names)
  expect_identical(rownames(V2), expected_names)
  expect_identical(colnames(V2), expected_names)
  expect_identical(rownames(s$coef_stage2), expected_names)

  reported_categories <- sub("^.*:", "", names(b2))

  expect_setequal(
    unique(reported_categories),
    expected_categories
  )
})
