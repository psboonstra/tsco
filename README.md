
# tsco

<!-- badges: start -->
<!-- badges: end -->

This R package fits two-stage conditional odds (tsco) models as described
in Boonstra, et al. (2026). Backend selection is automatic by default:
individual-level proportional-odds stages use `rms::orm()`, while multinomial
and grouped-count stages use `VGAM::vglm()`. Set `po_engine = "orm"` or
`po_engine = "vglm"` to explicitly request a proportional-odds engine; requesting
ORM for a grouped-count PO stage produces an error. A `tsco` object packages the
two fitted stage models together.

## Installation

You can install the development version of tsco from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("psboonstra/tsco")
```

## Example

``` r
library(tsco)
library(VGAM) # for wine data
data(wine)

# Proportional odds for Y | Y < bitter4
# Multinomial for {Y < bitter4, bitter4, bitter5}
tsco_ex <-
  tsco(
    formula = cbind(bitter1, bitter2, bitter3, bitter4, bitter5) ~ temp + contact,
    data = wine,
    cutoff_level = "bitter4",
    stage1 = "po",
    stage2 = "multinomial"
  )
 
#  
summary(tsco_ex)

predict(tsco_ex)
```

