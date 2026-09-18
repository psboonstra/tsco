#' Extract log-likelihood from a TSCO model
#'
#' Extracts the constant-free TSCO log-likelihood, computed on the
#' individual-level scale. For grouped/count responses, this excludes
#' multinomial combinatorial constants so that grouped and ungrouped
#' representations of the same data have comparable log-likelihoods.
#'
#' The `nobs` attribute, and therefore [BIC()], uses the frequency-weighted
#' sample size (`n_weighted`), so a weighted fit and a fit to the equivalent
#' row-replicated data give the same BIC. Without weights this is the number
#' of represented observations.
#'
#' @param object An object of class `"tsco"`.
#' @param ... Ignored.
#'
#' @return An object of class `"logLik"`.
#'
#' @method logLik tsco
#' @export
logLik.tsco <- function(object, ...) {
  val <- as.numeric(object$logLik)

  if (length(val) != 1L || !is.finite(val)) {
    stop("`object` does not contain a valid TSCO log-likelihood.", call. = FALSE)
  }

  df <- .tsco_df(object$fit_stage1) + .tsco_df(object$fit_stage2)

  nobs <- if (!is.null(object$n_weighted)) {
    object$n_weighted
  } else if (!is.null(object$n_obs)) {
    object$n_obs
  } else {
    object$n
  }

  structure(
    val,
    df = df,
    nobs = nobs,
    class = "logLik"
  )
}
