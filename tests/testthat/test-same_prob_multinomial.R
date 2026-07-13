test_that("tsco model produces same probabilities as standard multinomial in edge cases", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("dplyr")

  set.seed(100)

  n <- 500

  probs_x0 <- runif(5)
  probs_x1 <- runif(5)
  x <- rbinom(n, 1, 0.5)
  y <-
    (1 - x) * sample(5, size = n, prob = probs_x0, replace = TRUE) +
    x * sample(5, size = n, prob = probs_x1, replace = TRUE)

  dat <- data.frame(x = x, y = ordered(y, levels = c("1", "2", "3", "4", "5")))
  dat_grouped <- dat |>
    dplyr::group_by(x) |>
    dplyr::summarize(
      y1 = sum(y == "1"),
      y2 = sum(y == "2"),
      y3 = sum(y == "3"),
      y4 = sum(y == "4"),
      y5 = sum(y == "5"),
      .groups = "drop"
    )

  fit1 <- tsco(
    cbind(y1, y2, y3, y4, y5) ~ x,
    data = dat_grouped,
    cutoff_level = "y5",
    stage1 = "multinomial",
    stage2 = "po"
  )

  predict_fit1 <-
    predict(fit1, newdata = data.frame(x = c(0, 1)))

  fit2 <- tsco(
    cbind(y1, y2, y3, y4, y5) ~ x,
    data = dat_grouped,
    cutoff_level = "y5",
    stage1 = "multinomial",
    stage2 = "multinomial"
  )

  predict_fit2 <-
    predict(fit2, newdata = data.frame(x = c(0, 1)))


  fit3 <- VGAM::vglm(
    cbind(y1, y2, y3, y4, y5) ~ x,
    data = dat_grouped,
    family = VGAM::multinomial(),
  )

  predict_fit3 <-
    VGAM::predictvglm(fit3, newdata = data.frame(x = c(0, 1)), type = "response")

  expect_equal(predict_fit1, predict_fit2, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(predict_fit1, predict_fit3, tolerance = 1e-6, ignore_attr = TRUE)
})
