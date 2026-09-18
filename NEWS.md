# tsco 0.1.0

# tsco 0.0.0.9000

- `weights` of length `nrow(data)` are now passed through `model.frame()` as
  its `weights` argument, the contract `lm()` and `glm()` use, so they stay
  aligned with the retained rows under any `na.action`, including a user
  function that records omitted rows by name rather than position. The
  previous alignment through the `"na.action"` attribute assumed positions.
- **Breaking:** `tsco()` now errors early, naming the column, when a predictor
  in the model frame contains a missing or non-finite value, or the response
  is missing. This replaces opaque backend failures under
  `na.action = na.pass` and for transformations that are not finite on the
  data, such as `log(0)`.
- `predict()` treats non-finite predictor values (for example `log(x)` at
  `x = 0` or `x = Inf`) like missing ones: the row is returned as all-`NA`
  with a warning instead of a degenerate 0/1 probability row.
  Matrix-valued basis columns such as `poly()` are screened cell-wise.
- `predict(type = "stage")` now includes the `unfitted` and `incomplete` row
  indicators in its list.

- New `nobs()` method returning the frequency-weighted sample size
  `n_weighted`, matching the `nobs` attribute of `logLik()` and the sample
  size BIC uses. The documentation had referred to `nobs()` without the method
  existing.
- **Bug fix:** `predict()` with a missing value in a numeric predictor of
  `newdata` failed for the whole call ("could not be normalized" under
  `rms::orm()`, "unexpected number of rows" under `VGAM::vglm()`). Rows with
  any missing predictor are now returned as all-`NA`, with a warning and an
  `incomplete` attribute, and the remaining rows are predicted normally. A
  missing factor value is now reported the same way rather than as an
  unsupported factor level.
- **Breaking:** `data` is required and must be a data frame. The undocumented
  `parent.frame()` fallback is removed; it was untested and failed with
  `weights`.
- The `weights` documentation now states that weights are frequency weights
  throughout, including in standard errors and tests, and that survey (design)
  weights are out of scope.
- Comments and internal documentation no longer cite manuscript section or
  equation numbers, which may change; variants are named (`PO|C|MR`, etc.).

- Weighted `rms::orm()` stage fits no longer emit rms's note that weights are
  ignored by `validate()` and `bootcov()`. `tsco()` calls neither; the caveat
  stays in `?tsco` under `weights`. Only that one message is silenced.
- **Bug fix:** `weights` of length `nrow(data)` were aligned with the model
  frame by parsing `rownames()` as integers. Data frames with character row
  names errored, and numeric-looking row names that were not row positions
  selected the wrong weights silently. (The alignment now goes through
  `model.frame()`; see above.) Weights that are positive only on rows removed
  for missing data are now a clear error.
- **Breaking:** an `offset()` term in `formula` is now an error. Offsets were
  documented as unsupported but were still rebuilt into the stage formulas,
  where they failed inside `rms::orm()`.
- The estimate-only coefficient table returned when a stage's covariance
  cannot be extracted now carries the same reported names as `coef()`,
  `vcov()` and the full table, instead of the backend's raw names.

- **Breaking:** a formula without an intercept (`y ~ x - 1`, `y ~ 0 + x`) is
  now an error. The thresholds and multinomial intercepts are the baseline
  category probabilities and cannot be removed; previously `rms::orm()`
  ignored the `- 1` but switched to per-level slope names such as
  `temp=warm`, while `VGAM::vglm()` fitted a different multinomial model. The
  error message points to `relevel()` for choosing a factor's reference level.
- `predict(unsupported_stage1 = "uniform")` now splits the lower-partition
  mass over the lower categories observed in the fitting data. A category
  declared in `levels` but never observed receives probability zero, as it
  does under `"empirical"` and everywhere else in the package.

- **Weights are frequency weights, and BIC now says so.** The sample size used
  by BIC and by `nobs(logLik(fit))` is the weighted count `sum(weights)` (for
  grouped data, the weighted sum of the count totals), stored as `n_weighted`
  with the stage-1 analogue `n_stage1_weighted`, so a weighted fit and a fit
  to the equivalent row-replicated data now give the same BIC. Previously the
  log-likelihood was weighted but BIC used the unweighted row count. `n_obs`
  is unchanged. `print()` and `summary()` show the weighted counts only when
  weights were supplied.
- **Bug fix:** a row with case weight exactly zero is now removed before
  fitting. Such rows made `VGAM::vglm()` fail on a 0/0 in its log-likelihood,
  and would otherwise have counted as observed when deciding which predictor
  levels, outcome levels and stage-1 factor levels have support. `n_obs`,
  `n_weighted` and every support set now reflect positive-weight rows only;
  `n_zero_weight` records how many rows were removed. `predict()` without
  `newdata` still returns every input row; a row whose factor level occurred
  only in zero-weight observations is returned as `NA` with a warning, since
  no stage has information about it.
- Stage-2 fitting data now drop unused factor levels, as stage 1 already did,
  so `rms::orm()` no longer reports a zero coefficient with a singularity
  warning for a level that no observation has.
- The separation diagnostic now applies only to slope coefficients and also
  flags a finite estimate whose standard error is not finite. Threshold rows
  are excluded because `rms::vcov.orm()` may legitimately omit their
  covariance entries; that is not separation.
- Internal helpers no longer generate help pages (`@noRd`).

- **Bug fix:** the joint likelihood-ratio test for a single-predictor model
  fitted to a grouped-count (`cbind(...)`) response was badly
  anti-conservative. The full model's log-likelihood was taken from
  `VGAM::vglm()`, which includes the multinomial combinatorial constants for
  grouped data, while the closed-form intercept-only fit omits them; on the
  wine data the statistic was 162.2 where the individual-level fit gives 24.2.
  Multi-term LRTs were unaffected because the constants cancel between two
  refits. All LRTs now score both full and reduced models on the package's
  constant-free scale, so grouped and individual representations of the same
  data give the same test. Fitted objects gain `counts_stage1` and
  `counts_stage2`.
- **Bug fix:** `predict()` with raw `newdata` was silently wrong for formulas
  containing data-dependent basis terms such as `poly(x, 2)`, `splines::ns(x)`
  or `scale(x)`, because the basis was rebuilt from `newdata` alone. The stored
  terms now carry the fitted basis constants (`predvars`).
- **Bug fix:** `predict(unsupported_stage1 = "empirical")` fell back to a
  uniform split, with a warning, for a grouped-count fit whose `levels`
  declared an unobserved lower category, and ignored case weights in every
  case. It now uses the weighted stage-1 category totals and gives an
  unobserved level probability zero, consistent with the rest of the package.

- **Bug fix:** a level declared in `levels` but with no observations no longer
  breaks the fit. `rms::orm()` failed with an opaque
  `'names' attribute [5] must be the same length as the vector [4]`, raised
  inside `orm.fit()`; `VGAM::vglm()` silently dropped the level and returned a
  narrower probability matrix, which surfaced as
  `Predicted probability matrix has 4 columns, but 5 levels were expected`.
  Both stages, individual-level and grouped-count responses, and interior as
  well as trailing levels were affected. Such a level has a maximum likelihood
  fitted probability of exactly zero and a threshold on the boundary of the
  parameter space, so `tsco()` now drops it from the fit with a warning naming
  it, and `predict()` restores it as an all-zero column. Estimates and tests
  are conditional on the observed outcome levels; because the boundary sits in
  a threshold common to the full and reduced models, joint tests keep their
  usual degrees of freedom and chi-squared reference. A stage left with fewer
  than two observed levels is a clear error rather than a backend failure.
  Fitted objects gain `stage1_observed_levels`, `stage2_observed_levels` and
  `unobserved_levels`.
- Reported threshold names now follow the *observed* cutpoints. With an
  interior level unobserved the thresholds are `Y<=0`, `Y<=1`, `Y<=3`, not the
  first three declared levels.

- **Breaking:** reported coefficient names no longer depend on the fitting
  engine. Proportional-odds thresholds are named for the level they bound on
  the reported lower-tail scale (`Y<=2`, and `Y<C` for the stage-2 collapsed
  category) instead of `rms::orm()`'s upper-tail `y>=2` labels, which
  contradicted the sign-corrected value, or `VGAM::vglm()`'s positional
  `(Intercept):1`. PO slopes take design-matrix column names under both
  engines. `coef()`, `vcov()` and `summary()` are affected.
- `summary()` now reports per-stage convergence and separation diagnostics, and
  `print()` shows them when a stage fails to converge or a coefficient has the
  extreme-estimate-and-extreme-standard-error signature of complete or
  quasi-complete separation.
- Joint likelihood-ratio tests now work for single-predictor models such as
  `y ~ treatment`. The intercept-only reduced fit is computed in closed form
  from the observed category distribution, which such a model attains exactly,
  rather than by refitting a covariate-free model that not every backend
  supports.
- `weights` now works under the default `po_engine = "auto"`. Individual-level
  PO stages fall back to `VGAM::vglm()` when the installed `rms` does not
  accept case weights in `orm()`; `po_engine = "orm"` errors in that case
  rather than failing inside the fitter.
- Inline formula transformations such as `log(x)`, `poly(x, 2)` and
  `rms::rcs(x)` are now supported. (Offsets are not: `offset()` terms are
  carried into the stage formulas but fail inside `rms::orm()`.) Stage right-hand sides
  are rebuilt from model-frame column names, and `newdata` given to `predict()`
  may contain either the raw variables or the evaluated columns.
- `predict()` matches stage probability columns by name across both backends
  instead of silently trusting column order, and warns if it ever has to fall
  back to order.
- Printed coefficient tables now state the reported scale, and the component
  fit table notes that its per-stage BICs use per-stage sample sizes and so do
  not sum to the combined BIC.
- Fitted objects are substantially smaller: the stage fits no longer record a
  copy of the analysis data, the fitter function body, or the VGAM family
  object in their calls.
- Add a vignette on selecting a contextually motivated TsCO cutpoint and comparing conventional PO/MR and TsCO model structures.
- `logLik()` and `summary()` now use backend-reported effective model rank for degrees of freedom, with finite coefficient count as a fallback.
- `tsco()` adds `po_engine = "orm"` and `po_engine = "vglm"` to explicitly select the proportional-odds backend; the default `"auto"` uses `rms::orm()` for individual-level stages and `VGAM::vglm()` for grouped-count stages.
