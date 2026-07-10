test_that("sign correction is internally consistent between coef.tsco() and vcov.tsco()", {
  skip_if_not_installed("VGAM")

  # This does not test *which* sign convention is correct (see
  # test-sign_recovery_simulation.R and test-multinomial_sign_recovery.R for
  # that); it tests that whatever correction coef.tsco() applies is applied
  # identically -- and correctly, as a similarity transform D %*% V %*% D,
  # not a naive negation -- to vcov.tsco(). A mismatch here would mean
  # standard errors/CIs built from vcov() no longer correspond to the
  # coefficients coef() reports.

  data(pneumo, package = "VGAM")
  pneumo <- transform(pneumo, let = log(exposure.time))

  fit <- tsco(
    cbind(normal, mild, severe) ~ let,
    data = pneumo,
    levels = c("mild", "normal", "severe"),
    cutoff_level = "severe",
    stage1 = "po",
    stage2 = "multinomial",
    po.reverse = FALSE,
    warn_degenerate = FALSE
  )

  b1_raw <- stats::coef(fit$fit_stage1)
  b2_raw <- stats::coef(fit$fit_stage2)

  V1_raw <- as.matrix(stats::vcov(fit$fit_stage1))
  V2_raw <- as.matrix(stats::vcov(fit$fit_stage2))

  b1 <- coef(fit, stage = "stage1")
  b2 <- coef(fit, stage = "stage2")

  V1 <- vcov(fit, stage = "stage1")
  V2 <- vcov(fit, stage = "stage2")

  s1 <- unname(b1) / unname(b1_raw)
  s2 <- unname(b2) / unname(b2_raw)

  # Every ratio should be exactly +1 or -1 (a sign flip, not a rescaling).
  expect_true(all(abs(abs(s1) - 1) < 1e-8))
  expect_true(all(abs(abs(s2) - 1) < 1e-8))

  # Reconstruct the expected transformed covariance directly via outer(s, s)
  # and confirm it matches what vcov.tsco() actually returned -- i.e. that
  # vcov.tsco() used the same sign vector as coef.tsco(), not an
  # independently (and possibly inconsistently) derived one.
  V1_expected <- V1_raw * outer(s1, s1)
  V2_expected <- V2_raw * outer(s2, s2)

  expect_equal(unname(V1), unname(V1_expected), tolerance = 1e-10)
  expect_equal(unname(V2), unname(V2_expected), tolerance = 1e-10)

  # Sanity check the mechanics on a case where we know the answer by
  # construction: diagonal elements of V (variances) are invariant to any
  # +-1 sign flip, since (-1)^2 = 1. If this fails, the transform is
  # fundamentally broken (e.g. someone accidentally applied D %*% V, not
  # D %*% V %*% D).
  expect_equal(diag(V1), diag(V1_raw), tolerance = 1e-10)
  expect_equal(diag(V2), diag(V2_raw), tolerance = 1e-10)
})
