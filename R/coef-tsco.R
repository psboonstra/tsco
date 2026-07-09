#' Extract coefficients from a TSCO model
#'
#' @param object An object of class `"tsco"`.
#' @param stage Which stage to extract: `"both"`, `"stage1"`, or `"stage2"`.
#' @param ... Ignored.
#'
#' @return A named numeric vector of coefficients. If `stage = "both"`,
#'   coefficient names are prefixed by `"stage1:"` and `"stage2:"`.
#'
#' @method coef tsco
#' @export
coef.tsco <- function(object, stage = c("both", "stage1", "stage2"), ...) {
  stage <- match.arg(stage)

  b1 <- stats::coef(object$fit_stage1)
  b2 <- stats::coef(object$fit_stage2)

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
