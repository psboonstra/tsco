# A level declared in `levels` but absent from the data has a maximum
# likelihood fitted probability of exactly zero and a threshold on the boundary
# of the parameter space. Neither backend represents that: rms::orm() fails with
# an opaque `'names' attribute` error and VGAM::vglm() silently returns a
# narrower probability matrix. tsco() drops such levels from the fit, warns, and
# restores them at zero mass in predict().

make_sparse_data <- function(observed, n = 300, seed = 2) {
  set.seed(seed)
  data.frame(
    A = stats::rbinom(n, 1, 0.5),
    W = stats::rnorm(n),
    Y = sample(observed, n, replace = TRUE)
  )
}

all_levels <- as.character(0:5)

fit_sparse <- function(dat, cutoff_level = "5", ...) {
  tsco(
    Y ~ A + W,
    data = dat,
    levels = all_levels,
    cutoff_level = cutoff_level,
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE,
    ...
  )
}

test_that("an unobserved level is carried through at zero mass, either engine", {
  dat <- make_sparse_data(c("0", "1", "2", "3", "5"))

  for (engine in c("orm", "vglm")) {
    expect_warning(
      fit <- fit_sparse(dat, po_engine = engine),
      "level\\(s\\) 4 are declared in `levels` but have no observations"
    )

    expect_identical(fit$stage1_observed_levels, c("0", "1", "2", "3"))
    expect_identical(fit$unobserved_levels, "4")

    p <- suppressWarnings(predict(fit))

    expect_identical(colnames(p), all_levels)
    expect_equal(unname(p[, "4"]), rep(0, nrow(dat)))
    expect_equal(unname(rowSums(p)), rep(1, nrow(dat)), tolerance = 1e-12)
  }
})

test_that("thresholds are named for the observed cutpoints, not declared positions", {
  # The hazard case: with an interior level missing, VGAM's thresholds are the
  # cutpoints 0|1, 1|3 and 3|4. Naming them from positions 1..3 of the declared
  # set would label the 1|3 threshold `Y<=2` -- a fit that looks fine and is
  # silently mislabelled.
  dat <- make_sparse_data(c("0", "1", "3", "4", "5"))

  expect_warning(fit <- fit_sparse(dat), "level\\(s\\) 2 are declared")

  expect_identical(fit$stage1_observed_levels, c("0", "1", "3", "4"))

  expect_identical(
    names(coef(fit, stage = "stage1")),
    c("Y<=0", "Y<=1", "Y<=3", "A", "W")
  )

  expect_identical(
    rownames(vcov(fit, stage = "stage1")),
    names(coef(fit, stage = "stage1"))
  )

  p <- suppressWarnings(predict(fit))
  expect_equal(unname(p[, "2"]), rep(0, nrow(dat)))
})

test_that("the fit equals one that never declared the missing level", {
  dat <- make_sparse_data(c("0", "1", "2", "3", "5"))

  expect_warning(full <- fit_sparse(dat), "have no observations")

  reduced <- tsco(
    Y ~ A + W,
    data = dat,
    levels = setdiff(all_levels, "4"),
    cutoff_level = "5",
    stage1 = "po",
    stage2 = "po",
    warn_degenerate = FALSE
  )

  expect_equal(unname(coef(full)), unname(coef(reduced)), tolerance = 1e-8)
  expect_equal(
    as.numeric(logLik(full)), as.numeric(logLik(reduced)),
    tolerance = 1e-8
  )

  # The only difference is the extra all-zero column carried by the full fit.
  p_full <- suppressWarnings(predict(full))
  p_red <- predict(reduced)

  expect_identical(colnames(p_full), all_levels)
  expect_identical(colnames(p_red), setdiff(all_levels, "4"))
  strip <- function(m) matrix(as.numeric(m), nrow = nrow(m))

  expect_equal(
    strip(p_full[, colnames(p_red)]), strip(p_red),
    tolerance = 1e-8
  )
})

test_that("an unobserved level in stage 2 is handled the same way", {
  dat <- make_sparse_data(c("0", "1", "2", "3", "4"))

  expect_warning(
    fit <- fit_sparse(dat, cutoff_level = "4"),
    "Stage 2 .*level\\(s\\) 5 are declared"
  )

  expect_identical(fit$unobserved_levels, "5")

  p <- suppressWarnings(predict(fit))
  expect_equal(unname(p[, "5"]), rep(0, nrow(dat)))
  expect_equal(unname(rowSums(p)), rep(1, nrow(dat)), tolerance = 1e-12)
})

test_that("a grouped response with an all-zero column is handled the same way", {
  grouped <- data.frame(
    A = c(0, 1),
    `0` = c(20, 30), `1` = c(10, 5), `2` = c(8, 6),
    `3` = c(5, 9), `4` = c(0, 0), `5` = c(40, 20),
    check.names = FALSE
  )

  # Capture every warning, then assert on the package's own message. VGAM
  # warns about convergence at a half-step on this sparse table; backend
  # warnings are tolerated, but nothing else from the package may appear.
  warnings_seen <- character()
  fit <- withCallingHandlers(
    tsco(
      cbind(`0`, `1`, `2`, `3`, `4`, `5`) ~ A,
      data = grouped,
      levels = all_levels,
      cutoff_level = "5",
      stage1 = "po",
      stage2 = "po",
      warn_degenerate = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  is_ours <- grepl("declared in `levels` but have no observations", warnings_seen, fixed = TRUE)
  expect_equal(sum(is_ours), 1L)
  expect_match(warnings_seen[is_ours], "4", fixed = TRUE)

  others <- warnings_seen[!is_ours]
  backend <- grepl("half-step|convergence|iteration|checkwz", others)
  expect_length(others[!backend], 0L)

  expect_identical(fit$unobserved_levels, "4")

  p <- suppressWarnings(predict(fit))
  expect_identical(colnames(p), all_levels)
  expect_equal(unname(p[, "4"]), rep(0, nrow(grouped)))
  expect_equal(unname(rowSums(p)), rep(1, nrow(grouped)), tolerance = 1e-12)
})

test_that("fewer than two observed levels in a stage is a clear error", {
  dat <- make_sparse_data(c("0", "5"))

  expect_error(
    fit_sparse(dat),
    "has 1 outcome level\\(s\\) with observations"
  )

  expect_error(fit_sparse(dat), "At least two are needed")
})

test_that("nothing changes when every declared level is observed", {
  dat <- make_sparse_data(all_levels)

  expect_silent(fit <- fit_sparse(dat, po_engine = "vglm"))

  expect_identical(fit$stage1_observed_levels, fit$lower_levels)
  expect_identical(fit$stage2_observed_levels, fit$collapsed_levels)
  expect_identical(fit$unobserved_levels, character())

  p <- predict(fit)
  expect_identical(colnames(p), all_levels)
  expect_true(all(colSums(p) > 0))
})

test_that("the joint test keeps its degrees of freedom when a level is dropped", {
  # The boundary sits in a threshold, which appears identically in the full and
  # reduced models -- dropping a covariate does not change which outcome levels
  # are observed. Conditional on the observed support the tested coefficients
  # are interior, so the usual chi-square reference still applies.
  dat <- make_sparse_data(c("0", "1", "2", "3", "5"))

  expect_warning(fit <- fit_sparse(dat), "have no observations")

  lrt <- suppressWarnings(summary(fit, joint_test = "LRT")$joint_tests)

  expect_identical(rownames(lrt), c("A", "W"))
  expect_true(all(lrt$df == 2))
  expect_true(all(is.finite(lrt$Chisq)))
})
