#' Summarize a TSCO model
#'
#' Produces a concise statistical summary of both fitted TSCO stages.
#'
#' @param object An object of class `"tsco"`.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.tsco"`.
#'
#' @method summary tsco
#' @export
summary.tsco <- function(object, ...) {
  coef_stage1 <- .tsco_coef_table(object$fit_stage1, object$stage1)
  coef_stage2 <- .tsco_coef_table(object$fit_stage2, object$stage2)

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
    model = c(.tsco_model_label(object$stage1),
              .tsco_model_label(object$stage2)),
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
    stage1_label = .tsco_model_label(object$stage1),
    stage2_label = .tsco_model_label(object$stage2),
    po.reverse = object$po.reverse,
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
    coef_stage1 = coef_stage1,
    coef_stage2 = coef_stage2
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
    " using ", x$stage1_label, "\n",
    sep = ""
  )
  cat(
    "  Stage 2: collapsed outcome {Y < ", x$cutoff_level, ", ",
    paste(x$upper_levels, collapse = ", "),
    "} using ", x$stage2_label, "\n",
    sep = ""
  )

  cat("\nSample size:\n")
  if (isTRUE(x$grouped)) {
    cat("  Groups: ", x$n_groups, "\n", sep = "")
    cat("  Observations represented: ", x$n_obs, "\n", sep = "")
    cat("  Groups contributing to stage 1: ", x$n_stage1_groups, "\n", sep = "")
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

  cat("\nCombined fit statistics:\n")
  total_print <- x$total_fit
  for (jj in names(total_print)) {
    if (is.numeric(total_print[[jj]])) {
      total_print[[jj]] <- signif(total_print[[jj]], digits)
    }
  }
  print(total_print, row.names = FALSE, right = FALSE)

  cat("\nStage 1 coefficients: ")
  cat("Y | Y < ", x$cutoff_level, "\n", sep = "")

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

  invisible(x)
}
