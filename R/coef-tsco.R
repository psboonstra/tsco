#' Extract coefficients from a TSCO model
#'
#' Multinomial coefficients are reported in the manuscript's parameterization,
#' `Pr(Y = k | X) ~ exp(alpha_k - x'beta_k)`. For PO stages, covariate
#' coefficients use the manuscript's sign convention, while thresholds use the
#' common lower-tail convention
#' `Pr(Y <= j | X) = expit(gamma_j + x'beta)` across fitting engines. Relative
#' to manuscript cutpoints defined by
#' `Pr(Y >= k | X) = expit(alpha_k - x'beta)`, `gamma_j = -alpha_(j + 1)`.
#' Raw stage-fit coefficients are sign-corrected as needed; see
#' `.tsco_sign_vector()` for the derivation.
#'
#' Reported names do not depend on the fitting engine. Proportional-odds
#' thresholds are named `Y<=<level>` for the level they bound, matching the
#' lower-tail scale on which their values are reported, and slopes take the
#' design-matrix column names that `VGAM::vglm()` uses. In particular the
#' upper-tail threshold labels that `rms::orm()` produces (`y>=2`) are not
#' used, because they contradict the reported value once the sign correction
#' has been applied.
#'
#' @param object An object of class `"tsco"`.
#' @param stage Which stage to extract: `"both"`, `"stage1"`, or `"stage2"`.
#' @param ... Ignored.
#'
#' @return A named numeric vector on the reported coefficient scale described
#'   above. If `stage = "both"`, coefficient names are prefixed by `"stage1:"`
#'   and `"stage2:"`.
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
    engine = object$stage1_engine
  )

  s2 <- .tsco_sign_vector(
    object$fit_stage2,
    kind = object$stage2,
    engine = object$stage2_engine
  )

  b1 <- b1_raw * s1
  b2 <- b2_raw * s2

  names(b1) <- .tsco_apply_reported_names(
    names(b1_raw),
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

  names(b2) <- .tsco_apply_reported_names(
    names(b2_raw),
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
