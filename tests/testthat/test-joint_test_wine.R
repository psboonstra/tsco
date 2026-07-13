test_that("summary joint_test argument works", {
  dat <- make_wine_test_data()

  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  s_lrt <- summary(fit, joint_test = "LRT")
  s_wald <- summary(fit, joint_test = "Wald")
  s_none <- summary(fit, joint_test = "none")

  expect_true(nrow(s_lrt$joint_tests) > 0L)
  expect_true(nrow(s_wald$joint_tests) > 0L)
  expect_equal(nrow(s_none$joint_tests), 0L)

  expect_equal(s_lrt$joint_test_type, "LRT")
  expect_equal(s_wald$joint_test_type, "Wald")
  expect_equal(s_none$joint_test_type, "none")

  expect_error(summary(fit, joint_test = "banana"))
})
