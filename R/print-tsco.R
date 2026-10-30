#' Print a TSCO model
#'
#' @param x An object of class `"tsco"`.
#' @param ... Ignored.
#'
#' @method print tsco
#' @export
print.tsco <- function(x, ...) {
  cat("Two-Stage Conditional Odds Model\n")
  cat("--------------------------------\n")

  cat("Formula: ", .tsco_deparse(x$formula), "\n", sep = "")

  if (!is.null(x$grouped) && isTRUE(x$grouped)) {
    cat("Data type: grouped/count response\n")
  } else {
    cat("Data type: individual-level response\n")
  }

  cat("Outcome levels: ", paste(x$levels, collapse = " < "), "\n", sep = "")
  cat("Cutoff level: ", x$cutoff_level, "\n", sep = "")

  if (isTRUE(x$stage1_degenerate)) {
    cat(
      "Stage 1: Y | Y < ", x$cutoff_level, " is degenerate (",
      if (length(x$lower_levels) == 0L) "no categories" else "a single category",
      " below the cutoff); not fitted\n",
      sep = ""
    )
  } else {
    cat(
      "Stage 1: Y | Y < ", x$cutoff_level,
      " using ", .tsco_model_label(x$stage1), " [", x$stage1_engine, "]\n",
      sep = ""
    )
  }

  cat(
    "Stage 2: ",
    .tsco_stage2_outcome_label(x$cutoff_level, x$upper_levels, x$lower_levels),
    " using ", .tsco_model_label(x$stage2), " [", x$stage2_engine, "]\n",
    sep = ""
  )

  if (!is.null(x$grouped) && isTRUE(x$grouped)) {
    cat("Groups: ", x$n_groups, "\n", sep = "")
    cat("Observations represented: ", x$n_obs, "\n", sep = "")
    cat("Groups contributing to stage 1: ", x$n_stage1_groups, "\n", sep = "")
    cat("Observations contributing to stage 1: ", x$n_stage1_obs, "\n", sep = "")
  } else {
    cat("N: ", x$n_obs, "\n", sep = "")
    cat("N contributing to stage 1: ", x$n_stage1_obs, "\n", sep = "")
  }

  if (isTRUE(x$weighted)) {
    cat("Sum of weights: ", signif(x$n_weighted, 7), "\n", sep = "")
    cat(
      "Sum of weights contributing to stage 1: ",
      signif(x$n_stage1_weighted, 7), "\n",
      sep = ""
    )
  }

  cat("\nUse summary() for coefficient tables.\n")

  invisible(x)
}

