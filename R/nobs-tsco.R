#' Number of observations in a TSCO model
#'
#' Returns the frequency-weighted number of represented observations,
#' `n_weighted`: the sum of the case weights for an individual-level response,
#' or the weighted sum of the count totals for a grouped-count response. With
#' no weights this is the number of represented observations, `n_obs`. It is
#' the sample size that [BIC()] uses, and it matches the `nobs` attribute of
#' [logLik.tsco()].
#'
#' @param object An object of class `"tsco"`.
#' @param ... Ignored.
#'
#' @return A single number.
#'
#' @method nobs tsco
#' @importFrom stats nobs
#' @export
nobs.tsco <- function(object, ...) {
  if (!is.null(object$n_weighted)) {
    return(object$n_weighted)
  }

  if (!is.null(object$n_obs)) {
    return(object$n_obs)
  }

  object$n
}
