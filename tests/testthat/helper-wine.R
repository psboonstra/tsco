make_wine_test_data <- function() {
  skip_if_not_installed("ordinal")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")

  data(wine, package = "ordinal")

  wine <- as.data.frame(wine)

  # Ensure response is ordered and levels are character.
  wine$rating <- ordered(
    as.character(wine$rating),
    levels = levels(wine$rating)
  )

  lev <- levels(wine$rating)

  wine_grouped <- wine |>
    dplyr::count(temp, contact, rating, name = "n") |>
    tidyr::pivot_wider(
      names_from = rating,
      values_from = n,
      values_fill = 0
    ) |>
    as.data.frame()

  # Formula like cbind(`1`, `2`, `3`, `4`, `5`) ~ temp + contact
  grouped_formula <- stats::as.formula(
    paste0(
      "cbind(",
      paste(sprintf("`%s`", lev), collapse = ", "),
      ") ~ temp + contact"
    )
  )

  list(
    individual = wine,
    grouped = wine_grouped,
    levels = lev,
    cutoff_level = lev[4],
    individual_formula = rating ~ temp + contact,
    grouped_formula = grouped_formula
  )
}

fit_wine_individual <- function(stage1 = "po", stage2 = "po") {
  dat <- make_wine_test_data()

  tsco(
    dat$individual_formula,
    data = dat$individual,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = stage1,
    stage2 = stage2,
    warn_degenerate = FALSE
  )
}

fit_wine_grouped <- function(stage1 = "po", stage2 = "po") {
  dat <- make_wine_test_data()

  tsco(
    dat$grouped_formula,
    data = dat$grouped,
    levels = dat$levels,
    cutoff_level = dat$cutoff_level,
    stage1 = stage1,
    stage2 = stage2,
    warn_degenerate = FALSE
  )
}
