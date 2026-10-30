# A cutoff at the lowest or second outcome level leaves stage 1 with fewer than
# two categories. Its likelihood factor is then exactly 1, so the TSCO model is
# the stage-2 model on all K levels: a plain PO or multinomial fit reached
# through the same interface.

make_degenerate_data <- function(n = 300, seed = 3) {
  set.seed(seed)
  d <- data.frame(arm = factor(sample(c("control", "treatment"), n, TRUE)), x = stats::rnorm(n))
  lp <- 0.7 * (d$arm == "treatment") + 0.5 * d$x
  d$y <- ordered(cut(lp + stats::rlogis(n), c(-Inf, -0.7, 0.5, 1.5, Inf),
                     labels = c("a", "b", "c", "d")))
  d
}

to_grouped <- function(d) {
  tab <- table(d$arm, d$y)
  out <- data.frame(arm = factor(rownames(tab)), as.data.frame.matrix(tab), check.names = FALSE)
  rownames(out) <- NULL
  out
}


test_that("a cutoff at the lowest or second level fits stage 2 only, with a warning", {
  d <- make_degenerate_data()

  for (cl in c("a", "b")) {
    expect_warning(
      fit <- tsco(y ~ arm + x, data = d, cutoff_level = cl, stage1 = "po", stage2 = "po"),
      "is degenerate"
    )

    expect_true(fit$stage1_degenerate)
    expect_null(fit$fit_stage1)
    expect_null(fit$data_stage1)
    expect_null(fit$counts_stage1)
    expect_null(fit$totals_stage1)
    expect_identical(fit$stage1_engine, NA_character_)
    expect_null(fit$stage1_formula)
    expect_equal(fit$logLik_stage1, 0)
    expect_length(fit$stage1_dropped_terms, 0L)
    expect_length(coef(fit, stage = "stage1"), 0L)
    expect_equal(dim(vcov(fit, stage = "stage1")), c(0L, 0L))

    expect_no_warning(
      tsco(y ~ arm + x, data = d, cutoff_level = cl, stage1 = "po", stage2 = "po",
           warn_degenerate = FALSE)
    )
  }

  # The lowest level has no collapsed category; the second does.
  fa <- suppressWarnings(tsco(y ~ arm + x, data = d, cutoff_level = "a", warn_degenerate = FALSE))
  fb <- suppressWarnings(tsco(y ~ arm + x, data = d, cutoff_level = "b", warn_degenerate = FALSE))
  expect_null(fa$lower_collapsed_label)
  expect_identical(fa$collapsed_levels, c("a", "b", "c", "d"))
  expect_identical(fa$lower_levels, character())
  expect_identical(fb$lower_levels, "a")
  expect_identical(fb$collapsed_levels, c(fb$lower_collapsed_label, "b", "c", "d"))
})


test_that("the degenerate model equals a plain PO fit on all levels", {
  d <- make_degenerate_data()
  ref <- VGAM::vglm(y ~ arm + x, data = d,
                    family = VGAM::cumulative(link = "logitlink", parallel = TRUE))
  ref_ll <- sum(log(VGAM::predictvglm(ref, type = "response")[cbind(seq_len(nrow(d)), as.integer(d$y))]))

  for (cl in c("a", "b")) {
    fit <- suppressWarnings(tsco(y ~ arm + x, data = d, cutoff_level = cl,
                                 stage1 = "po", stage2 = "po", warn_degenerate = FALSE))

    expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-5)
    expect_equal(as.numeric(logLik(fit)), ref_ll, tolerance = 1e-6)
    expect_equal(attr(logLik(fit), "df"), 5)
    expect_equal(nobs(fit), nrow(d))
    expect_equal(BIC(fit), -2 * ref_ll + log(nrow(d)) * 5, tolerance = 1e-6)

    p <- predict(fit, newdata = d[1:5, c("arm", "x")])
    p_ref <- VGAM::predictvglm(ref, newdata = d[1:5, ], type = "response")
    expect_equal(unclass(p)[, ], unname(p_ref), tolerance = 1e-5, ignore_attr = TRUE)
    expect_false(any(attr(p, "unsupported_stage1")))
    expect_equal(unname(rowSums(p)), rep(1, 5), tolerance = 1e-8)

    cls <- predict(fit, newdata = d[1:5, c("arm", "x")], type = "class")
    expect_identical(as.character(cls), colnames(p_ref)[max.col(p_ref)])

    st <- predict(fit, newdata = d[1:5, c("arm", "x")], type = "stage")
    expect_equal(dim(st$stage1_conditional), c(5L, length(fit$lower_levels)))
    expect_equal(st$prob, p, ignore_attr = TRUE)
    if (cl == "a") expect_equal(st$lower_mass, rep(0, 5))
    if (cl == "b") expect_equal(st$lower_mass, unname(p_ref[, "a"]), tolerance = 1e-5)
  }

  # The reported names differ only in how the first threshold is labeled.
  fa <- suppressWarnings(tsco(y ~ arm + x, data = d, cutoff_level = "a", warn_degenerate = FALSE, stage2 = "po"))
  fb <- suppressWarnings(tsco(y ~ arm + x, data = d, cutoff_level = "b", warn_degenerate = FALSE, stage2 = "po"))
  expect_identical(names(coef(fa))[1:2], c("stage2:Y<=a", "stage2:Y<=b"))
  expect_identical(names(coef(fb))[1:2], c("stage2:Y<b", "stage2:Y<=b"))
})


test_that("joint tests for a degenerate stage 1 are the stage-2 tests", {
  d <- make_degenerate_data()
  fit <- suppressWarnings(tsco(y ~ arm + x, data = d, cutoff_level = "b",
                               stage1 = "po", stage2 = "po", warn_degenerate = FALSE))

  full <- VGAM::vglm(y ~ arm + x, data = d, family = VGAM::cumulative(parallel = TRUE))
  no_arm <- VGAM::vglm(y ~ x, data = d, family = VGAM::cumulative(parallel = TRUE))
  no_x <- VGAM::vglm(y ~ arm, data = d, family = VGAM::cumulative(parallel = TRUE))
  ll <- function(f) sum(log(VGAM::predictvglm(f, type = "response")[cbind(seq_len(nrow(d)), as.integer(d$y))]))

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  expect_identical(rownames(lrt), c("arm", "x"))
  expect_equal(lrt$df, c(1, 1))
  expect_equal(lrt["arm", "Chisq"], 2 * (ll(full) - ll(no_arm)), tolerance = 1e-5)
  expect_equal(lrt["x", "Chisq"], 2 * (ll(full) - ll(no_x)), tolerance = 1e-5)

  wald <- summary(fit, joint_test = "Wald")$joint_tests
  expect_identical(rownames(wald), c("arm", "x"))
  expect_equal(wald$df, c(1, 1))
  b <- coef(fit); V <- vcov(fit)
  expect_equal(wald["arm", "Chisq"], unname(b["stage2:armtreatment"]^2 / V["stage2:armtreatment", "stage2:armtreatment"]), tolerance = 1e-8)
})


test_that("degenerate stage 1 works for multinomial stage 2 and grouped data", {
  d <- make_degenerate_data()
  g <- to_grouped(d)

  fit <- suppressWarnings(tsco(cbind(a, b, c, d) ~ arm, data = g, cutoff_level = "a",
                               stage2 = "multinomial", warn_degenerate = FALSE))
  expect_true(fit$stage1_degenerate)
  expect_identical(fit$stage2_engine, "vglm")

  full <- suppressWarnings(VGAM::vglm(cbind(a, b, c, d) ~ arm, data = g, family = VGAM::multinomial(refLevel = 1)))
  null <- suppressWarnings(VGAM::vglm(cbind(a, b, c, d) ~ 1, data = g, family = VGAM::multinomial(refLevel = 1)))
  counts <- as.matrix(g[, c("a", "b", "c", "d")])
  ll <- function(f) sum(counts * log(VGAM::predictvglm(f, type = "response")))

  # Multinomial slopes are reported on the manuscript's alpha - x'beta scale,
  # i.e. negated relative to VGAM; intercepts are unchanged.
  ref_coef <- coef(full)
  is_int <- grepl("(Intercept)", names(ref_coef), fixed = TRUE)
  ref_coef[!is_int] <- -ref_coef[!is_int]
  expect_equal(unname(coef(fit)), unname(ref_coef), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(fit)), ll(full), tolerance = 1e-6)
  expect_equal(attr(logLik(fit), "df"), 6)

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  expect_equal(lrt$df, 3)
  expect_equal(lrt$Chisq, 2 * (ll(full) - ll(null)), tolerance = 1e-6)

  s <- summary(fit, joint_test = "none")
  expect_identical(s$component_fit$model[1], "degenerate")
  expect_true(is.na(s$component_fit$BIC[1]))
  expect_output(print(s), "Degenerate stage; no parameters")
  expect_output(print(s), "Reference category: Y = a")
  expect_output(print(fit), "Stage 2: \\{a, b, c, d\\}")
  expect_output(print(fit), "degenerate \\(no categories below the cutoff\\)")

  # A grouped fit with the cutoff at the second level, PO stage 2.
  fb <- suppressWarnings(tsco(cbind(a, b, c, d) ~ arm, data = g, cutoff_level = "b",
                              stage2 = "po", warn_degenerate = FALSE))
  ref <- VGAM::vglm(cbind(a, b, c, d) ~ arm, data = g, family = VGAM::cumulative(parallel = TRUE))
  expect_equal(unname(coef(fb)), unname(coef(ref)), tolerance = 1e-6)
  expect_output(print(fb), "degenerate \\(a single category below the cutoff\\)")
})


test_that("the degenerate warning names the stage-2 model type and how to get PO", {
  d <- make_degenerate_data()

  expect_warning(
    tsco(y ~ arm + x, data = d, cutoff_level = "a"),
    "multinomial model on all 4 outcome levels \\(`stage2 = \"multinomial\"`; for a proportional-odds model set `stage2 = \"po\"`\\)"
  )
  expect_warning(
    tsco(y ~ arm + x, data = d, cutoff_level = "a", stage2 = "po"),
    "proportional odds model on all 4 outcome levels \\(`stage2 = \"po\"`\\)"
  )
})


test_that("a non-degenerate cutoff is unaffected", {
  d <- make_degenerate_data()
  fit <- tsco(y ~ arm + x, data = d, cutoff_level = "c", stage1 = "po", stage2 = "po")
  expect_false(fit$stage1_degenerate)
  expect_false(is.null(fit$fit_stage1))
  expect_true(any(grepl("^stage1:", names(coef(fit)))))
})
