
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

This is a basic example which shows you how to solve a common problem:

``` r
library(tsco)
wine_fmla <- as.formula(cbind(bitter1, bitter2, bitter3, bitter4, bitter5) ~ temp + contact)

# Multinomial model using VGAM
vglm(wine_fmla, family = multinomial, data = wine)

# Proportional odds model using VGAM
vglm(wine_fmla, family = cumulative(link = "logitlink", parallel = T, reverse = F), data = wine)

#
tsco(
  wine_fmla,
  data = wine,
  cutoff_level = "bitter4",
  stage1 = "po",
  stage2 = "po"
)
 
```

