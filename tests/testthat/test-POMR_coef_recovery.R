test_that("tsco recovers known true multinomial-stage coefficients (simple, 1-param case)", {
  skip_if_not_installed("VGAM")

  set.seed(20260710)

  n <- 5e3
  x <- stats::rnorm(n)

  true_beta1 <- 0.8    # stage 1 (conditional PO): L1 vs L2
  true_beta_C <- -1.5  # stage 2 (marginal MR), single non-reference category

  alpha1 <- 0.4
  alpha_C <- 0.6

  # --- Stage 2: two categories only, "<C" (reference) vs "C" --------------
  # Pr(Y2 = "<C") ~ 1 ; Pr(Y2 = "C") ~ exp(alpha_C - x * true_beta_C)
  odds_C <- exp(alpha_C - x * true_beta_C)
  p2_C <- odds_C / (1 + odds_C)

  u2 <- stats::runif(n)
  y2 <- ifelse(u2 <= p2_C, "C", "<C")

  # --- Stage 1: only meaningful for subjects with y2 == "<C" --------------
  p1_ge2 <- stats::plogis(alpha1 - x * true_beta1) # Pr(Y1 >= "L2")
  u1 <- stats::runif(n)
  y1 <- ifelse(u1 <= (1 - p1_ge2), "L1", "L2")

  y <- ifelse(y2 == "<C", y1, y2)

  dat <- data.frame(x = x, y = ordered(y, levels = c("L1", "L2", "C")))

  fit <- tsco(
    y ~ x,
    data = dat,
    levels = c("L1", "L2", "C"),
    cutoff_level = "C",
    stage1 = "po",
    stage2 = "multinomial",
    po.reverse = FALSE,
    warn_degenerate = FALSE
  )

  b2 <- coef(fit, stage = "stage2")
  x_coef2 <- b2[grepl("^x$|:x$|^x:", names(b2)) | names(b2) == "x"]

  expect_length(x_coef2, 1L)

  expect_equal(unname(x_coef2), true_beta_C, tolerance = 0.15)
  expect_true(sign(unname(x_coef2)) == sign(true_beta_C))
})


test_that("tsco recovers known true multinomial-stage coefficients (2-param case, order-assumption documented)", {
  skip_if_not_installed("VGAM")

  set.seed(20260710)

  n <- 5e3
  x <- stats::rnorm(n)

  true_beta1 <- 0.8
  true_beta_C  <- 1.0   # deliberately different magnitude/sign from U2
  true_beta_U2 <- -0.5  # so a naive sign-flip or swap is distinguishable

  alpha1 <- 0.4
  alpha_C  <- 0.6
  alpha_U2 <- -0.2

  # --- Stage 2: three categories, "<C" (reference), "C", "U2" -------------
  odds_C  <- exp(alpha_C  - x * true_beta_C)
  odds_U2 <- exp(alpha_U2 - x * true_beta_U2)
  denom <- 1 + odds_C + odds_U2

  p2_lt_C <- 1 / denom
  p2_C    <- odds_C / denom
  # p2_U2 <- odds_U2 / denom, implied

  u2 <- stats::runif(n)
  y2 <- ifelse(
    u2 <= p2_lt_C, "<C",
    ifelse(u2 <= p2_lt_C + p2_C, "C", "U2")
  )

  # --- Stage 1: only meaningful for subjects with y2 == "<C" --------------
  p1_ge2 <- stats::plogis(alpha1 - x * true_beta1)
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
    stage2 = "multinomial",
    po.reverse = FALSE,
    warn_degenerate = FALSE
  )

  b2 <- coef(fit, stage = "stage2")
  x_coef2 <- b2[grepl("^x$|:x$|^x:", names(b2)) | names(b2) == "x"]

  expect_length(x_coef2, 2L)

  # First diagnostic: print names so a test-output review can immediately
  # confirm/deny the ordering assumption stated above.
  # (testthat will show this on failure via the surrounding expectations,
  # but an explicit label makes it unambiguous.)
  x_coef2_named <- x_coef2
  if (length(x_coef2_named) == 2L) {
    names(x_coef2_named) <- names(x_coef2)
  }

  expect_equal(
    unname(x_coef2[1]), true_beta_C,
    tolerance = 0.15,
    label = paste0("coef named '", names(x_coef2)[1], "' assumed to be 'C'")
  )
  expect_equal(
    unname(x_coef2[2]), true_beta_U2,
    tolerance = 0.15,
    label = paste0("coef named '", names(x_coef2)[2], "' assumed to be 'U2'")
  )

  expect_true(sign(unname(x_coef2[1])) == sign(true_beta_C))
  expect_true(sign(unname(x_coef2[2])) == sign(true_beta_U2))
})
