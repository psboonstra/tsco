#' Variance-covariance matrix for a TSCO model
#'
#' Reported on the same sign-corrected scale as `coef.tsco()`. If a stage's
#' coefficients required a sign flip to match the manuscript's
#' parameterization (see `.tsco_sign_vector()`), the covariance matrix is
#' transformed as `V_new = D V_old D`, where `D` is the diagonal +-1 sign
#' matrix -- not simply negated -- so that cross terms between flipped and
#' unflipped coefficients (e.g. an intercept and a slope) remain correct.
#'
#' @param object An object of class `"tsco"`.
#' @param stage Which stage to extract: `"both"`, `"stage1"`, or `"stage2"`.
#' @param ... Ignored.
#'
#' @return A variance-covariance matrix on the TSCO/manuscript coefficient
#'   scale. If `stage = "both"`, the returned matrix is block diagonal.
#'
#' @method vcov tsco
#' @export
vcov.tsco <- function(object, stage = c("both", "stage1", "stage2"), ...) {
  stage <- match.arg(stage)

  V1_raw <- as.matrix(stats::vcov(object$fit_stage1))
  V2_raw <- as.matrix(stats::vcov(object$fit_stage2))

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

  # Align sign vectors to covariance matrix names when possible.
  s1 <- .tsco_align_sign_vector(s1, V1_raw)
  s2 <- .tsco_align_sign_vector(s2, V2_raw)

  V1 <- V1_raw * outer(s1, s1)
  V2 <- V2_raw * outer(s2, s2)

  if (stage == "stage1") {
    return(V1)
  }

  if (stage == "stage2") {
    return(V2)
  }

  p1 <- nrow(V1)
  p2 <- nrow(V2)

  i1 <- seq_len(p1)
  i2 <- p1 + seq_len(p2)

  V <- matrix(0, nrow = p1 + p2, ncol = p1 + p2)

  V[i1, i1] <- V1
  V[i2, i2] <- V2

  rn1 <- rownames(V1)
  rn2 <- rownames(V2)

  if (is.null(rn1)) {
    rn1 <- names(stats::coef(object$fit_stage1))
  }

  if (is.null(rn2)) {
    rn2 <- names(stats::coef(object$fit_stage2))
  }

  rn <- c(
    paste0("stage1:", rn1),
    paste0("stage2:", rn2)
  )

  dimnames(V) <- list(rn, rn)

  V
}
