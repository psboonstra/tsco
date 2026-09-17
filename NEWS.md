# tsco 0.0.0.9000

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
  `rms::rcs(x)` are now supported, along with offsets. Stage right-hand sides
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
