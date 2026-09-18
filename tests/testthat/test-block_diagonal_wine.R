test_that("combined vcov is block diagonal", {
  dat <- make_wine_test_data()

  # This is the manuscript's headline PO|C|MR variant: conditional PO
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

  V <- vcov(fit)
  V1 <- vcov(fit, stage = "stage1")
  V2 <- vcov(fit, stage = "stage2")

  p1 <- nrow(V1)
  p2 <- nrow(V2)

  expect_equal(
    V[seq_len(p1), seq_len(p1)],
    V1,
    ignore_attr = TRUE
  )

  expect_equal(
    V[p1 + seq_len(p2), p1 + seq_len(p2)],
    V2,
    ignore_attr = TRUE
  )

  expect_equal(
    V[seq_len(p1), p1 + seq_len(p2)],
    matrix(0, p1, p2),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )

  expect_equal(
    V[p1 + seq_len(p2), seq_len(p1)],
    matrix(0, p2, p1),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})
