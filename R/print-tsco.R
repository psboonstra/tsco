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

  cat(
    "Stage 1: Y | Y < ", x$cutoff_level,
    " using ", .tsco_model_label(x$stage1), " [", x$stage1_engine, "]\n",
    sep = ""
  )

  cat(
    "Stage 2: {Y < ", x$cutoff_level, ", ",
    paste(x$upper_levels, collapse = ", "),
    "} using ", .tsco_model_label(x$stage2), " [", x$stage2_engine, "]\n",
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

  cat("\nUse summary() for coefficient tables.\n")

  invisible(x)
}

