# Variables whose names need backticks (`test arm`), including inside inline
# transformations (log(`my x`)). The stage formulas are rebuilt from
# model-frame column names, so quoting has to be exact, and term matching for
# the joint tests has to compare canonical variable sets rather than labels.

make_odd_names_data <- function(n = 240, seed = 11) {
  set.seed(seed)
  d <- data.frame(
    check.names = FALSE,
    `test arm` = factor(sample(c("control", "treatment"), n, TRUE)),
    `my x` = exp(stats::runif(n, 0, 2))
  )
  lp <- 0.8 * (d$`test arm` == "treatment") + 0.5 * log(d$`my x`)
  d$y <- ordered(
    cut(lp + stats::rlogis(n), c(-Inf, -0.5, 0.8, 2, Inf), labels = c("a", "b", "c", "d"))
  )
  d
}


test_that("a non-syntactic factor name fits, tests and predicts (the reported reproducer)", {
  skip_if_not_installed("tibble")

  grouped <- tibble::tibble(
    `test arm` = factor(c("control", "treatment")),
    a = c(18, 14), b = c(14, 12), c = c(10, 18), d = c(3, 4)
  )

  fit <- tsco(cbind(a, b, c, d) ~ `test arm`, data = grouped, cutoff_level = "d",
              warn_degenerate = FALSE)

  expect_identical(attr(stats::terms(fit$stage1_formula), "term.labels"), "`test arm`")
  expect_true("stage1:`test arm`treatment" %in% names(coef(fit)))

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  wald <- summary(fit, joint_test = "Wald")$joint_tests
  expect_identical(rownames(lrt), "`test arm`")
  expect_identical(rownames(wald), "`test arm`")
  expect_equal(lrt$df, 2)   # 1 slope in stage 1 + 1 in the binary MR stage 2
  expect_true(is.finite(lrt$Chisq))

  nd <- data.frame(check.names = FALSE, `test arm` = factor(c("control", "treatment")))
  p <- predict(fit, newdata = nd)
  expect_equal(dim(p), c(2L, 4L))
  expect_equal(unname(rowSums(p)), c(1, 1), tolerance = 1e-8)
  expect_equal(nrow(predict(fit)), 2L)
})


test_that("odd names give the same model as syntactic names, including transformations", {
  d <- make_odd_names_data()
  d2 <- d
  names(d2) <- c("arm", "x", "y")

  fit <- tsco(y ~ `test arm` * log(`my x`), data = d, cutoff_level = "c",
              stage1 = "po", stage2 = "po", warn_degenerate = FALSE)
  ref <- tsco(y ~ arm * log(x), data = d2, cutoff_level = "c",
              stage1 = "po", stage2 = "po", po_engine = "vglm", warn_degenerate = FALSE)

  # The factor name is non-syntactic, so PO stages go to VGAM.
  expect_identical(fit$stage1_engine, "vglm")
  expect_identical(fit$stage2_engine, "vglm")

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-6)

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  lrt_ref <- summary(ref, joint_test = "LRT")$joint_tests
  expect_identical(rownames(lrt), c("`test arm`", "log(`my x`)", "`test arm`:log(`my x`)"))
  expect_equal(lrt$df, lrt_ref$df)
  expect_equal(lrt$df, c(4, 4, 2))
  expect_equal(lrt$Chisq, lrt_ref$Chisq, tolerance = 1e-6)

  wald <- summary(fit, joint_test = "Wald")$joint_tests
  wald_ref <- summary(ref, joint_test = "Wald")$joint_tests
  expect_identical(rownames(wald), rownames(lrt))
  expect_equal(wald$df, c(4, 4, 2))
  expect_equal(wald$Chisq, wald_ref$Chisq, tolerance = 1e-6)

  p <- predict(fit, newdata = d[1:5, c("test arm", "my x")])
  p_ref <- predict(ref, newdata = d2[1:5, c("arm", "x")])
  expect_equal(unclass(p)[, ], unclass(p_ref)[, ], tolerance = 1e-6, ignore_attr = TRUE)
})


test_that("explicit po_engine = 'orm' with a non-syntactic factor name is a clear error", {
  d <- make_odd_names_data()
  expect_error(
    tsco(y ~ `test arm`, data = d, cutoff_level = "c", stage1 = "po", stage2 = "po",
         po_engine = "orm", warn_degenerate = FALSE),
    "non-syntactic names \\(`test arm`\\)"
  )

  # A non-syntactic *numeric* name is fine for orm.
  fit <- tsco(y ~ log(`my x`), data = d, cutoff_level = "c", stage1 = "po", stage2 = "po",
              po_engine = "orm", warn_degenerate = FALSE)
  expect_identical(fit$stage1_engine, "orm")
})


test_that("the two-arm constant-predictor drop works with a non-syntactic name", {
  d <- make_odd_names_data()
  d$y[d$`test arm` == "treatment" & d$y %in% c("a", "b")] <- "c"

  expect_warning(
    fit <- tsco(y ~ `test arm`, data = d, cutoff_level = "c", stage1 = "po",
                stage2 = "po", warn_degenerate = FALSE),
    "`test arm` take a single value"
  )
  expect_identical(fit$stage1_dropped_terms, "`test arm`")
  expect_length(attr(stats::terms(fit$stage1_formula), "term.labels"), 0L)

  lrt <- summary(fit, joint_test = "LRT")$joint_tests
  expect_identical(rownames(lrt), "`test arm`")
  expect_equal(lrt$df, 1)
  expect_true(is.finite(lrt$Chisq))
})


test_that("canonical term matching sees through quoting", {
  user <- stats::terms(y ~ `test arm` * log(`my x`))
  stage <- stats::as.formula("z ~ `test arm` + `log(\\`my x\\`)` + `test arm`:`log(\\`my x\\`)`")

  expect_identical(
    .tsco_term_variables(user)[["log(`my x`)"]],
    .tsco_term_variables(stats::terms(stage))[["`log(\\`my x\\`)`"]]
  )
  expect_setequal(
    .tsco_terms_containing("log(`my x`)", stage, user),
    c("`log(\\`my x\\`)`", "`test arm`:`log(\\`my x\\`)`")
  )
  expect_setequal(
    .tsco_terms_containing("`test arm`", stage, user),
    c("`test arm`", "`test arm`:`log(\\`my x\\`)`")
  )
  expect_identical(.tsco_terms_containing("nope", stage, user), character())

  expect_identical(.tsco_quote_name(c("arm", "test arm", "log(`my x`)")),
                   c("arm", "`test arm`", "`log(\\`my x\\`)`"))
  expect_identical(.tsco_canonical_name(c("`test arm`", "`log(\\`my x\\`)`", "arm")),
                   c("test arm", "log(`my x`)", "arm"))
})
