test_that("tsco recovers known true coefficients (sign and magnitude) under po.reverse = TRUE", {
  skip_if_not_installed("VGAM")

  # -------------------------------------------------------------------------
  # Data-generating mechanism.
  #
  # We simulate directly from the TSCO factorization using the manuscript's
  # parameterization: Pr(Y >= k | X = x) = expit(alpha_k - x * beta).
  #
  # Four ordinal levels: "L1" < "L2" < "C" < "U2".
  # cutoff_level = "C", so:
  #   lower_levels     = c("L1", "L2")        (stage 1: conditional PO)
  #   collapsed_levels = c("<C", "C", "U2")   (stage 2: marginal PO)
  #
  # Stage 2 (marginal, always observed) determines whether a subject falls
  # below the cutoff or lands in "C"/"U2". Stage 1 (conditional, only defined
  # for subjects below the cutoff) then determines "L1" vs "L2".
  #
  # If tsco()'s po.reverse = TRUE default matches this parameterization, the
  # fitted x-coefficient in each stage should recover true_beta1 / true_beta2
  # in both sign and magnitude. A same-magnitude, flipped-sign result would
  # indicate VGAM's `reverse` argument does not mean what the documentation
  # in tsco.R claims.
  # -------------------------------------------------------------------------

  set.seed(20260709)

  n <- 50000
  x <- stats::rnorm(n)

  true_beta1 <- 0.8   # stage 1 (conditional PO on Y < C): L1 vs L2
  true_beta2 <- -1.2  # stage 2 (marginal PO): <C vs C vs U2

  alpha1 <- 0.4        # single cutpoint for a 2-level stage-1 outcome
  alpha2 <- c(1.0, -0.3) # two cutpoints for a 3-level stage-2 outcome

  # --- Stage 2 draw: assigns every subject to "<C", "C", or "U2" ----------
  p2_ge2 <- stats::plogis(alpha2[1] - x * true_beta2) # Pr(Y2 >= "C")
  p2_ge3 <- stats::plogis(alpha2[2] - x * true_beta2) # Pr(Y2 >= "U2")

  # Enforce proper ordering (should hold given alpha2 chosen above, but
  # guard against numerical edge cases at extreme x).
  p2_ge3 <- pmin(p2_ge3, p2_ge2)

  u2 <- stats::runif(n)
  y2 <- ifelse(u2 <= (1 - p2_ge2), "<C",
               ifelse(u2 <= (1 - p2_ge3), "C", "U2"))

  # --- Stage 1 draw: only meaningful for subjects with y2 == "<C" --------
  p1_ge2 <- stats::plogis(alpha1 - x * true_beta1) # Pr(Y1 >= "L2")
  u1 <- stats::runif(n)
  y1 <- ifelse(u1 <= (1 - p1_ge2), "L1", "L2")

  y <- ifelse(y2 == "<C", y1, y2)

  dat <- data.frame(x = x, y = ordered(y, levels = c("L1", "L2", "C", "U2")))

  fit <- tsco(
    y ~ x,
    data = dat,
    levels = c("L1", "L2", "C", "U2"),
    cutoff_level = "C",
    stage1 = "po",
    stage2 = "po",
    po.reverse = TRUE,
    warn_degenerate = FALSE
  )

  b1 <- coef(fit, stage = "stage1")
  b2 <- coef(fit, stage = "stage2")

  # Identify the x coefficient in each stage robustly (name is whatever VGAM
  # assigns; assume it is not one of the intercept/cutpoint terms).
  x_coef1 <- b1[grepl("^x$|:x$|^x:", names(b1)) | names(b1) == "x"]
  x_coef2 <- b2[grepl("^x$|:x$|^x:", names(b2)) | names(b2) == "x"]

  expect_length(x_coef1, 1L)
  expect_length(x_coef2, 1L)

  # With n = 50000, SEs on this coefficient should be small (~0.02-0.03);
  # tolerance is generous relative to that to avoid flakiness while still
  # catching a sign error (which would miss by ~2 * true_beta, not ~0.1).
  expect_equal(unname(x_coef1), true_beta1, tolerance = 0.15)
  expect_equal(unname(x_coef2), true_beta2, tolerance = 0.15)

  # Explicit sign checks, stated separately so a failure here is
  # unambiguous about *what* went wrong (sign vs. magnitude).
  expect_true(sign(unname(x_coef1)) == sign(true_beta1))
  expect_true(sign(unname(x_coef2)) == sign(true_beta2))
})
