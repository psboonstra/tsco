#' Summarize a TSCO model
#'
#' Produces a concise statistical summary of both fitted TSCO stages.
#'
#' Likelihood-ratio tests can be preferable to Wald tests when coefficient
#' estimates are unstable, but they do not resolve complete or quasi-complete
#' separation. If a stage's maximum-likelihood estimate is attained only at
#' infinity, both Wald inference and the ordinary chi-squared reference
#' distribution for likelihood-ratio tests may be unreliable. Inspect
#' convergence diagnostics and coefficient estimates before interpreting either
#' test.
#'
#' The returned object carries a `diagnostics` element reporting, per stage,
#' whether the backend converged and whether any coefficient shows the
#' magnitude-and-standard-error signature of complete or quasi-complete
#' separation. `print()` shows these only when there is something to report.
#' They are a warning that the numbers need scrutiny, not a correction.
#'
#' @param object An object of class `"tsco"`.
#' @param joint_test Character string specifying the type of joint test to
#'  perform across stages. Options are `"LRT"` for likelihood ratio test,
#'  `"Wald"` for Wald test, or `"none"` to skip joint tests. Default is `"LRT"`.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.tsco"`.
#'
#' @method summary tsco
#' @export
summary.tsco <- function(object, joint_test = c("LRT", "Wald", "none"), ...) {
  coef_stage1 <- .tsco_coef_table(
    object$fit_stage1,
    kind = object$stage1,
    engine = object$stage1_engine,
    reported_names = .tsco_stage_reported_names(
      object$fit_stage1,
      kind = object$stage1,
      engine = object$stage1_engine,
      response_levels = .tsco_stage_levels(object, 1L),
      formula = object$stage1_formula,
      data = object$data_stage1,
      collapsed_label = object$lower_collapsed_label,
      cutoff_level = object$cutoff_level
    )
  )

  coef_stage2 <- .tsco_coef_table(
    object$fit_stage2,
    kind = object$stage2,
    engine = object$stage2_engine,
    reported_names = .tsco_stage_reported_names(
      object$fit_stage2,
      kind = object$stage2,
      engine = object$stage2_engine,
      response_levels = .tsco_stage_levels(object, 2L),
      formula = object$stage2_formula,
      data = object$data_stage2,
      collapsed_label = object$lower_collapsed_label,
      cutoff_level = object$cutoff_level
    )
  )

  joint_test <- match.arg(joint_test)

  joint_tests <- switch(
    joint_test,
    LRT = .tsco_joint_lrt_tests(object),
    Wald = .tsco_joint_wald_tests(object),
    none = data.frame()
  )

  ll1 <- as.numeric(object$logLik_stage1)
  ll2 <- as.numeric(object$logLik_stage2)

  df1 <- .tsco_df(object$fit_stage1)
  df2 <- .tsco_df(object$fit_stage2)

  aic1 <- if (is.finite(ll1) && is.finite(df1)) -2 * ll1 + 2 * df1 else NA_real_
  aic2 <- if (is.finite(ll2) && is.finite(df2)) -2 * ll2 + 2 * df2 else NA_real_

  bic1 <- if (is.finite(ll1) && is.finite(df1)) {
    -2 * ll1 + log(object$n_stage1_obs) * df1
  } else {
    NA_real_
  }

  bic2 <- if (is.finite(ll2) && is.finite(df2)) {
    -2 * ll2 + log(object$n_obs) * df2
  } else {
    NA_real_
  }

  component_fit <- data.frame(
    stage = c("Stage 1", "Stage 2"),
    model = c(
      .tsco_model_label(object$stage1),
      .tsco_model_label(object$stage2)
    ),
    engine = c(object$stage1_engine, object$stage2_engine),
    outcome = c(
      paste0("Y | Y < ", object$cutoff_level),
      paste0("{Y < ", object$cutoff_level, ", ",
             paste(object$upper_levels, collapse = ", "), "}")
    ),
    n_obs = c(object$n_stage1_obs, object$n_obs),
    n_groups = c(object$n_stage1_groups, object$n_groups),
    df = c(df1, df2),
    logLik = c(ll1, ll2),
    AIC = c(aic1, aic2),
    BIC = c(bic1, bic2),
    stringsAsFactors = FALSE
  )

  total_logLik <- as.numeric(object$logLik)
  total_df <- if (is.finite(df1) && is.finite(df2)) df1 + df2 else NA_real_

  total_AIC <- if (is.finite(total_logLik) && is.finite(total_df)) {
    -2 * total_logLik + 2 * total_df
  } else {
    NA_real_
  }

  total_BIC <- if (is.finite(total_logLik) && is.finite(total_df)) {
    -2 * total_logLik + log(object$n_obs) * total_df
  } else {
    NA_real_
  }

  out <- list(
    call = object$call,
    formula = object$formula,
    grouped = isTRUE(object$grouped),
    levels = object$levels,
    cutoff_level = object$cutoff_level,
    cutoff_index = object$cutoff_index,
    lower_levels = object$lower_levels,
    upper_levels = object$upper_levels,
    collapsed_levels = object$collapsed_levels,
    lower_collapsed_label = object$lower_collapsed_label,
    stage1 = object$stage1,
    stage2 = object$stage2,
    stage1_engine = object$stage1_engine,
    stage2_engine = object$stage2_engine,
    stage1_label = .tsco_model_label(object$stage1),
    stage2_label = .tsco_model_label(object$stage2),
    n = object$n_obs,
    n_obs = object$n_obs,
    n_groups = object$n_groups,
    n_stage1 = object$n_stage1_obs,
    n_stage1_obs = object$n_stage1_obs,
    n_stage1_groups = object$n_stage1_groups,
    component_fit = component_fit,
    total_fit = data.frame(
      n_obs = object$n_obs,
      df = total_df,
      logLik = total_logLik,
      AIC = total_AIC,
      BIC = total_BIC
    ),
    joint_test_type = joint_test,
    joint_tests = joint_tests,
    coef_stage1 = coef_stage1,
    coef_stage2 = coef_stage2,
    diagnostics = list(
      stage1 = .tsco_stage_diagnostics(
        object$fit_stage1,
        engine = object$stage1_engine,
        coef_table = coef_stage1
      ),
      stage2 = .tsco_stage_diagnostics(
        object$fit_stage2,
        engine = object$stage2_engine,
        coef_table = coef_stage2
      )
    )
  )

  class(out) <- "summary.tsco"
  out
}


#' Print a TSCO summary
#'
#' @param x An object of class `"summary.tsco"`.
#' @param digits Number of digits to print.
#' @param signif.stars Logical; whether to print significance stars.
#' @param ... Additional arguments passed to print methods.
#'
#' @method print summary.tsco
#' @export
print.summary.tsco <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    signif.stars = getOption("show.signif.stars"),
    ...) {

  cat("Two-Stage Conditional Odds Model\n")
  cat("================================\n\n")

  cat("Formula: ", .tsco_deparse(x$formula), "\n", sep = "")

  if (isTRUE(x$grouped)) {
    cat("Data type: grouped/count response\n")
  } else {
    cat("Data type: individual-level response\n")
  }

  cat("Outcome levels: ", paste(x$levels, collapse = " < "), "\n", sep = "")
  cat("Cutoff level: ", x$cutoff_level, "\n", sep = "")
  cat("Lower partition: ", paste(x$lower_levels, collapse = ", "), "\n", sep = "")
  cat("Upper partition: ", paste(x$upper_levels, collapse = ", "), "\n", sep = "")

  cat("\nModel structure:\n")
  cat(
    "  Stage 1: Y | Y < ", x$cutoff_level,
    " using ", x$stage1_label, " [", x$stage1_engine, "]\n",
    sep = ""
  )
  cat(
    "  Stage 2: collapsed outcome {Y < ", x$cutoff_level, ", ",
    paste(x$upper_levels, collapse = ", "),
    "} using ", x$stage2_label, " [", x$stage2_engine, "]\n",
    sep = ""
  )

  cat("\nSample size:\n")
  if (isTRUE(x$grouped)) {
    cat("  Unique covariate patterns: ", x$n_groups, "\n", sep = "")
    cat("  Observations represented: ", x$n_obs, "\n", sep = "")
    cat("  Unique covariate patterns contributing to stage 1: ", x$n_stage1_groups, "\n", sep = "")
    cat("  Observations contributing to stage 1: ", x$n_stage1_obs, "\n", sep = "")
  } else {
    cat("  N: ", x$n_obs, "\n", sep = "")
    cat("  N contributing to stage 1: ", x$n_stage1_obs, "\n", sep = "")
  }

  cat("\nComponent fit statistics:\n")
  fit_print <- x$component_fit

  if (!isTRUE(x$grouped) && "n_groups" %in% names(fit_print)) {
    fit_print$n_groups <- NULL
  }

  numeric_cols <- c("n_obs", "n_groups", "df", "logLik", "AIC", "BIC")
  for (jj in intersect(numeric_cols, names(fit_print))) {
    fit_print[[jj]] <- signif(fit_print[[jj]], digits)
  }
  print(fit_print, row.names = FALSE, right = FALSE)
  cat(
    "  Note: each stage's BIC uses that stage's own sample size, so the two\n",
    "  component BICs do not sum to the combined BIC below.\n",
    sep = ""
  )

  cat("\nCombined fit statistics:\n")
  total_print <- x$total_fit
  for (jj in names(total_print)) {
    if (is.numeric(total_print[[jj]])) {
      total_print[[jj]] <- signif(total_print[[jj]], digits)
    }
  }
  print(total_print, row.names = FALSE, right = FALSE)

  if (x$joint_test_type != "none") {
    cat("\nJoint", x$joint_test_type, "tests across stages:\n")
    cat("  H0: all stage-specific coefficients for the term are zero\n")
    if (is.null(x$joint_tests) || nrow(x$joint_tests) == 0L) {
      cat("  No joint tests available.\n")
    } else {
      stats::printCoefmat(
        as.matrix(x$joint_tests),
        digits = digits,
        signif.stars = signif.stars,
        na.print = "NA",
        has.Pvalue = TRUE,
        P.values = TRUE,
        ...
      )
    }
  } else {
    cat("\nNo joint tests performed.\n")
  }

  cat("\nStage 1 coefficients: ")
  cat("Y | Y < ", x$cutoff_level, "\n", sep = "")
  cat(.tsco_scale_legend(x$stage1), sep = "")
  if (identical(x$stage1, "multinomial")) {
    cat(
      "Reference category: Y = ",
      x$lower_levels[1],
      "\n",
      sep = ""
    )
  }


  if (length(x$coef_stage1) == 0L) {
    cat("  No coefficient table available.\n")
  } else {
    stats::printCoefmat(
      x$coef_stage1,
      digits = digits,
      signif.stars = signif.stars,
      na.print = "NA",
      ...
    )
  }

  cat("\nStage 2 coefficients: ")
  cat("collapsed outcome {Y < ", x$cutoff_level, ", ",
      paste(x$upper_levels, collapse = ", "), "}\n", sep = "")
  cat(.tsco_scale_legend(x$stage2), sep = "")
  if (identical(x$stage2, "multinomial")) {
    cat(
      "Reference category: Y < ",
      x$cutoff_level,
      "\n",
      sep = ""
    )
  }

  if (length(x$coef_stage2) == 0L) {
    cat("  No coefficient table available.\n")
  } else {
    stats::printCoefmat(
      x$coef_stage2,
      digits = digits,
      signif.stars = signif.stars,
      na.print = "NA",
      ...
    )
  }

  .tsco_print_diagnostics(x$diagnostics)

  invisible(x)
}

#' One-line reminder of the reported coefficient scale for a stage.
#'
#' Both stage parameterizations subtract the linear predictor, following the
#' manuscript, so a positive coefficient shifts probability toward *lower*
#' outcome categories. That is the opposite of what most users expect from a
#' regression table, and the two stages label their intercepts differently, so
#' the convention is stated next to every printed table.
#'
#' @param kind "po" or "multinomial".
.tsco_scale_legend <- function(kind) {
  if (identical(kind, "multinomial")) {
    return(paste0(
      "Scale: Pr(Y = k | X) proportional to exp(alpha_k - x'beta_k); ",
      "intercepts are alpha_k.\n",
      "       A positive beta shifts probability toward lower categories.\n"
    ))
  }

  paste0(
    "Scale: Pr(Y >= k | X) = expit(alpha_k - x'beta); thresholds are printed ",
    "on the\n       lower-tail scale, Y<=j meaning gamma_j = -alpha_(j+1).\n",
    "       A positive beta shifts probability toward lower categories.\n"
  )
}

#' Print convergence and separation diagnostics, when there are any.
#'
#' @param diagnostics The `diagnostics` element of a `summary.tsco` object.
.tsco_print_diagnostics <- function(diagnostics) {
  if (is.null(diagnostics)) {
    return(invisible(NULL))
  }

  lines <- character()

  for (stage in names(diagnostics)) {
    d <- diagnostics[[stage]]

    if (is.null(d) || length(d$messages) == 0L) {
      next
    }

    label <- if (identical(stage, "stage1")) "Stage 1" else "Stage 2"
    lines <- c(lines, paste0("  ", label, ": ", d$messages))
  }

  if (length(lines) == 0L) {
    return(invisible(NULL))
  }

  cat("\nFit diagnostics:\n")
  cat(lines, sep = "\n")
  cat(
    "\n  Wald tests and the chi-squared calibration of likelihood-ratio tests",
    "\n  may be unreliable for the affected coefficients.\n"
  )

  invisible(NULL)
}
