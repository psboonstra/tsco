#' Variance-covariance matrix for a TSCO model
#'
#' Reported on the same sign-corrected scale as `coef.tsco()`: manuscript-scale
#' multinomial coefficients, and manuscript-sign covariate coefficients plus
#' common lower-tail thresholds for PO stages. If coefficients require a sign
#' flip (see `.tsco_sign_vector()`), the covariance matrix is transformed as
#' `V_new = D V_old D`, where `D` is the diagonal +-1 sign matrix -- not simply
#' negated -- so that cross terms between flipped and unflipped coefficients
#' remain correct.
#'
#' For `rms::orm()` stages, all threshold covariances are requested. If an
#' `orm` fit stores only a reduced covariance matrix and the complete matrix
#' cannot be recovered, unavailable threshold rows and columns are returned as
#' `NA`; slope covariance entries remain available when supplied by the fit.
#' Joint Wald tests exclude threshold rows.
#'
#' @param object An object of class `"tsco"`.
#' @param stage Which stage to extract: `"both"`, `"stage1"`, or `"stage2"`.
#' @param ... Ignored.
#'
#' @return A variance-covariance matrix on the reported coefficient scale
#'   described above. If `stage = "both"`, the returned matrix is block
#'   diagonal.
#'
#' @method vcov tsco
#' @export
vcov.tsco <- function(object, stage = c("both", "stage1", "stage2"), ...) {
  stage <- match.arg(stage)

  V1_raw <- .tsco_stage_vcov(object$fit_stage1)
  V2_raw <- .tsco_stage_vcov(object$fit_stage2)

  if (is.null(V1_raw) || is.null(V2_raw)) {
    stop("Could not extract a covariance matrix for each TSCO stage.", call. = FALSE)
  }

  s1 <- .tsco_sign_vector(
    object$fit_stage1,
    kind = object$stage1,
    engine = object$stage1_engine
  )

  s2 <- .tsco_sign_vector(
    object$fit_stage2,
    kind = object$stage2,
    engine = object$stage2_engine
  )

  # Align sign vectors to covariance matrix names when possible.
  s1 <- .tsco_align_sign_vector(s1, V1_raw)
  s2 <- .tsco_align_sign_vector(s2, V2_raw)

  V1 <- V1_raw * outer(s1, s1)
  V2 <- V2_raw * outer(s2, s2)

  rn1 <- .tsco_apply_reported_names(
    rownames(V1),
    .tsco_stage_reported_names(
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

  rn2 <- .tsco_apply_reported_names(
    rownames(V2),
    .tsco_stage_reported_names(
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

  dimnames(V1) <- list(rn1, rn1)
  dimnames(V2) <- list(rn2, rn2)

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
