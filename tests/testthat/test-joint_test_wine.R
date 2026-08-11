test_that("Wald tests combine all coefficients belonging to a formula term", {
  dat <- data.frame(
    grp = factor(rep(c("A", "B", "C"), each = 80)),
    z = factor(rep(rep(c("low", "high"), each = 40), 3)),
    y = ordered(
      rep(c("a", "b", "c", "d"), 60),
      levels = c("a", "b", "c", "d")
    )
  )

  for (engine in c("auto", "vglm")) {
    fit <- tsco(
      y ~ grp + z + grp:z,
      data = dat,
      cutoff_level = "c",
      stage1 = "po",
      stage2 = "po",
      po_engine = engine,
      warn_degenerate = FALSE
    )

    wald <- summary(fit, joint_test = "Wald")$joint_tests

    expect_identical(rownames(wald), c("grp", "z", "grp:z"))
    expect_equal(wald$df, c(4, 2, 4))
    expect_equal(sum(!is.finite(wald$Chisq)), 0L)
  }
})


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
