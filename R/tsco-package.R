#' tsco: Two-Stage Conditional Odds Models
#'
#' The `tsco` package fits two-stage conditional odds (TSCO) models for
#' ordinal outcomes. TSCO models partition an ordinal response at a chosen
#' cutoff level and fit two component models:
#'
#' * a conditional model for the lower partition, `Y | Y < C`;
#' * a marginal model for the collapsed outcome `{Y < C, C, ..., K - 1}`.
#'
#' The package currently uses `VGAM::vglm()` to fit the component models.
#' Each component may be specified as a proportional odds model or a
#' multinomial model.
#'
#' @section Main function:
#'
#' The primary user-facing function is [tsco()], which fits a TSCO model.
#'
#' @section Supported response formats:
#'
#' `tsco()` supports both individual-level ordinal outcomes, for example
#'
#' ```
#' y ~ x1 + x2
#' ```
#'
#' and grouped count outcomes, for example
#'
#' ```
#' cbind(y0, y1, y2) ~ x1 + x2
#' ```
#'
#' For grouped count responses, the optional `levels` argument can be used to
#' specify the ordinal ordering of the response categories.
#'
#' @section Model variants:
#'
#' The arguments `stage1` and `stage2` determine the two component model types.
#' Current options are:
#'
#' * `"po"`: proportional odds;
#' * `"multinomial"`: multinomial baseline-category logit.
#'
#' For example, `stage1 = "po"` and `stage2 = "multinomial"` fits a
#' `PO | C | MR` TSCO model.
#'
#' @section Package conventions:
#'
#' The package reports TSCO log-likelihoods on a constant-free individual-level
#' scale. This makes grouped and ungrouped representations of the same data
#' comparable.
#'
#' @seealso [tsco()], [summary.tsco()], [predict.tsco()]
#'
#' @docType package
#' @name tsco-package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
