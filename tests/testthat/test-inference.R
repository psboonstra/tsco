test_that("summary reports coefficient SEs and p-values", {
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

  b <- coef(fit)
  V <- vcov(fit)

  base_names <- .tsco_base_coef_name(names(b))
  idx <- which(base_names == "let")

  b_let <- b[idx]
  V_let <- V[names(b_let), names(b_let), drop = FALSE]

  stat_manual <- as.numeric(crossprod(b_let, solve(V_let, b_let)))

  s <- summary(fit)

  expect_true(all(
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in%
      colnames(s$coef_stage1)
  ))

  expect_true(all(
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in%
      colnames(s$coef_stage2)
  ))

  expect_true("let" %in% rownames(s$joint_tests))
  expect_equal(s$joint_tests["let", "df"], 2)
  expect_true(is.finite(s$joint_tests["let", "Chisq"]))
  expect_true(is.finite(s$joint_tests["let", "Pr(>Chisq)"]))

  expect_equal(
    s$joint_tests["let", "Chisq"],
    stat_manual,
    tolerance = 1e-8
  )

})
