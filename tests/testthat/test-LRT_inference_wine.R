test_that("LRT in summary matches manual reduced-model refit", {
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

  s <- summary(fit, joint_test = "LRT")

  expect_true(all(
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in%
      colnames(s$coef_stage1)
  ))

  expect_true(all(
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in%
      colnames(s$coef_stage2)
  ))

  expect_true(all(c("temp", "contact") %in% rownames(s$joint_tests)))
  expect_true(all(s$joint_tests[c("temp", "contact"), "df"] == 2))
  expect_true(all(is.finite(s$joint_tests[c("temp", "contact"), "Chisq"])))
  expect_true(all(is.finite(s$joint_tests[c("temp", "contact"), "Pr(>Chisq)"])))

  expect_true(all(is.finite(s$joint_tests$Chisq)))
  expect_true(all(s$joint_tests$Chisq >= -1e-8))
  expect_true(all(is.finite(s$joint_tests[["Pr(>Chisq)"]])))

  term <- attr(stats::terms(fit$stage1_formula), "term.labels")[2]

  f1_reduced <- .tsco_drop_term_formula(fit$stage1_formula, term)
  f2_reduced <- .tsco_drop_term_formula(fit$stage2_formula, term)

  fam1 <- .tsco_make_family(fit$stage1, po.reverse = fit$po.reverse)
  fam2 <- .tsco_make_family(fit$stage2, po.reverse = fit$po.reverse)

  fit1_reduced <- .tsco_fit_vglm_internal(
    formula = f1_reduced,
    family = fam1,
    data = fit$data_stage1,
    weights = fit$weights_stage1,
    vglm_args = fit$vglm_args
  )

  fit2_reduced <- .tsco_fit_vglm_internal(
    formula = f2_reduced,
    family = fam2,
    data = fit$data_stage2,
    weights = fit$weights_stage2,
    vglm_args = fit$vglm_args
  )

  ll_full <- .tsco_vglm_loglik(fit$fit_stage1) +
    .tsco_vglm_loglik(fit$fit_stage2)

  ll_reduced <- .tsco_vglm_loglik(fit1_reduced) +
    .tsco_vglm_loglik(fit2_reduced)

  lrt_manual <- 2 * (ll_full - ll_reduced)

  df_manual <- (
    .tsco_df(fit$fit_stage1) - .tsco_df(fit1_reduced)
  ) + (
    .tsco_df(fit$fit_stage2) - .tsco_df(fit2_reduced)
  )

  expect_equal(
    s$joint_tests[term, "Chisq"],
    lrt_manual,
    tolerance = 1e-8
  )

  expect_equal(
    s$joint_tests[term, "df"],
    df_manual
  )

})
