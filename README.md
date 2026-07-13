
# tsco

<!-- badges: start -->
<!-- badges: end -->

This is R package fits two-stage conditional odds (tsco) models as described
in Boonstra, et al. (2026). Much of the heavy lifting is done by VGAM, and the 
code follows the general semantic structure of `vglm`. A new class of model
`tsco` is created, which is essentially two `vglm` objects appropriately packaged
together. 

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

# Multinomial for Y | Y < bitter4
# Proportional odds for {Y < bitter 4, bitter 4, bitter 5}
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

