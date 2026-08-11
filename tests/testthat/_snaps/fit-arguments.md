# stage-specific fitter arguments are validated

    Code
      tsco(dat$individual_formula, data = dat$individual, levels = dat$levels,
      cutoff_level = dat$cutoff_level, stage1_args = TRUE)
    Condition
      Error:
      ! `stage1_args` must be a list.

---

    Code
      tsco(dat$individual_formula, data = dat$individual, levels = dat$levels,
      cutoff_level = dat$cutoff_level, stage1_args = list(TRUE))
    Condition
      Error:
      ! Every element of `stage1_args` must be named.

---

    Code
      tsco(dat$individual_formula, data = dat$individual, levels = dat$levels,
      cutoff_level = dat$cutoff_level, stage2_args = structure(list(TRUE, FALSE),
      names = c("x.arg", "x.arg")))
    Condition
      Error:
      ! Duplicated argument names in `stage2_args`: x.arg.

---

    Code
      tsco(dat$individual_formula, data = dat$individual, levels = dat$levels,
      cutoff_level = dat$cutoff_level, stage2_args = list(data = dat$individual,
      family = "logistic"))
    Condition
      Error:
      ! `stage2_args` cannot contain arguments controlled by `tsco()`: data, family.

