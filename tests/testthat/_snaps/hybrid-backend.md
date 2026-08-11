# grouped PO stages retain the VGAM backend

    Code
      tsco(cbind(a, b, c, d) ~ x, data = dat, cutoff_level = "c", stage1 = "po",
      stage2 = "po", po_engine = "orm", warn_degenerate = FALSE)
    Condition
      Error:
      ! `po_engine = "orm"` cannot be used with grouped-count PO stages; use `po_engine = "auto"` or `po_engine = "vglm"`.

