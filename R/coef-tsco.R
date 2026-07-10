#' Extract coefficients from a TSCO model
#'
#' Coefficients are reported in the manuscript's parameterization
#' (Pr(Y >= k | X) = expit(alpha_k - x'beta) for PO stages;
#' Pr(Y = k | X) ~ exp(alpha_k - x'beta_k) for multinomial stages), which
#' may require sign-correcting the raw VGAM fit's coefficients; see
#' `.tsco_sign_vector()` for the derivation.
#'
#' @param object An object of class `"tsco"`.
#' @param stage Which stage to extract: `"both"`, `"stage1"`, or `"stage2"`.
#' @param ... Ignored.
#'
#' @return A named numeric vector of coefficients on the TSCO/manuscript scale.
#'   If `stage = "both"`, coefficient names are prefixed by `"stage1:"` and
#'   `"stage2:"`.
#'
#' @method coef tsco
#' @export
coef.tsco <- function(object, stage = c("both", "stage1", "stage2"), ...) {
  stage <- match.arg(stage)

  b1_raw <- stats::coef(object$fit_stage1)
  b2_raw <- stats::coef(object$fit_stage2)

  s1 <- .tsco_sign_vector(
    object$fit_stage1,
    kind = object$stage1,
    po.reverse = object$po.reverse
  )

  s2 <- .tsco_sign_vector(
    object$fit_stage2,
    kind = object$stage2,
    po.reverse = object$po.reverse
  )

  b1 <- b1_raw * s1
  b2 <- b2_raw * s2

  if (stage == "stage1") {
    return(b1)
  }

  if (stage == "stage2") {
    return(b2)
  }

  names(b1) <- paste0("stage1:", names(b1))
  names(b2) <- paste0("stage2:", names(b2))

  c(b1, b2)
}
