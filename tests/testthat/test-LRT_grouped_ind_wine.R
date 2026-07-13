test_that("grouped and individual wine LRTs agree", {
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

  sg <- summary(fit_grouped, joint_test = "LRT")
  si <- summary(fit_ind, joint_test = "LRT")

  expect_equal(
    rownames(sg$joint_tests),
    rownames(si$joint_tests)
  )

  expect_equal(
    sg$joint_tests$df,
    si$joint_tests$df
  )

  expect_equal(
    sg$joint_tests$Chisq,
    si$joint_tests$Chisq,
    tolerance = 1e-7
  )

  expect_equal(
    sg$joint_tests[["Pr(>Chisq)"]],
    si$joint_tests[["Pr(>Chisq)"]],
    tolerance = 1e-7
  )
})
