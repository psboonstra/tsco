# tsco 0.0.0.9000

- Add a vignette on selecting a contextually motivated TsCO cutpoint and comparing conventional PO/MR and TsCO model structures.
- `logLik()` and `summary()` now use backend-reported effective model rank for degrees of freedom, with finite coefficient count as a fallback.
- `tsco()` adds `po_engine = "orm"` and `po_engine = "vglm"` to explicitly select the proportional-odds backend; the default `"auto"` uses `rms::orm()` for individual-level stages and `VGAM::vglm()` for grouped-count stages.
