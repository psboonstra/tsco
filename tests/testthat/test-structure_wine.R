test_that("PO|C|MR (stage1 = po, stage2 = multinomial) fits and is well-formed", {
  dat <- make_wine_test_data()

  # This is the manuscript's headline variant (Section 2.1.1): conditional PO
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

  expect_s3_class(fit, "tsco")
  expect_true(is.finite(fit$logLik))
  expect_true(is.finite(fit$logLik_stage1))
  expect_true(is.finite(fit$logLik_stage2))

  b <- coef(fit)
  V <- vcov(fit)

  expect_equal(length(b), nrow(V))
  expect_equal(names(b), rownames(V))

  # vcov should still be block diagonal regardless of stage2's family.
  p1 <- length(coef(fit, stage = "stage1"))
  p2 <- length(coef(fit, stage = "stage2"))
  expect_equal(
    V[seq_len(p1), p1 + seq_len(p2)],
    matrix(0, p1, p2),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )

  # print()/summary() should not error for a non-PO/PO combination.
  expect_output(print(fit))
  expect_output(print(summary(fit)))
})

test_that("MR|C|PO (stage1 = multinomial, stage2 = po) fits and is well-formed", {
  dat <- make_wine_test_data()

  # The mirror-image variant flagged in Section 2.1.3 as a natural extension:
  # conditional (saturated) multinomial below the cutoff, marginal PO at/above.
  fit <- tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = "multinomial",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_s3_class(fit, "tsco")
  expect_true(is.finite(fit$logLik))

  b <- coef(fit)
  V <- vcov(fit)

  expect_equal(length(b), nrow(V))
  expect_equal(names(b), rownames(V))

  expect_output(print(fit))
  expect_output(print(summary(fit)))
})

test_that("all four stage1/stage2 combinations produce parameter counts matching the manuscript's formulas", {
  dat <- make_wine_test_data()

  # p = 2 predictor (temp, contact, both binary). C (cutoff_index) = 4 here: lower_levels has 3
  # categories (1, 2, 3), upper_levels has 2 (4, 5), K = 5.
  p <- 2L
  K <- 5L
  cutoff_index <- 4L
  n_lower <- cutoff_index - 1L # = 3
  n_upper <- K - cutoff_index + 1L # = 2

  fit_combo <- function(stage1, stage2) {
    fit <- tsco(
      dat$grouped_formula,
      data = dat$grouped,
      levels = dat$levels,
      cutoff_level = dat$cutoff_level,
      stage1 = stage1,
      stage2 = stage2,
      warn_degenerate = FALSE
    )
  }

  # --- PO|C|PO: K - 1 + 2p total (manuscript, Section 2.1.2) --------------
  fit_po_po <- fit_combo("po", "po")
  expect_equal(
    length(coef(fit_po_po)),
    (K - 1L) + 2L * p
  )

  # --- PO|C|MR: (C - 1 + p) + (K - C)(p + 1) (Section 2.1.1) --------------
  # Here "C" in the manuscript's formula is cutoff_index.
  fit_po_mr <- fit_combo("po", "multinomial")
  expect_equal(
    length(coef(fit_po_mr)),
    (n_lower - 1L + p) + (n_upper) * (p + 1L)
  )

  # --- MR|C|PO: symmetric construction, derived the same way -------------
  # Stage 1 (multinomial, saturated over n_lower categories):
  #   (n_lower - 1)(p + 1) parameters.
  # Stage 2 (PO, over n_upper + 1 collapsed categories): n_upper + p params.
  fit_mr_po <- fit_combo("multinomial", "po")
  expect_equal(
    length(coef(fit_mr_po)),
    (n_lower - 1L) * (p + 1L) + (n_upper + p)
  )

  # --- MR|C|MR: symmetric construction, derived the same way -------------
  # Stage 1 (multinomial, saturated over n_lower categories):
  #   (n_lower - 1)(p + 1) parameters.
  # Stage 2 (multinomial, over n_upper + 1 collapsed categories): n_upper * (p + 1)
  fit_mr_mr <- fit_combo("multinomial", "multinomial")
  expect_equal(
    length(coef(fit_mr_mr)),
    (n_lower - 1L) * (p + 1L) + (n_upper * (p + 1L))
  )
})
