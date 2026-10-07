---
title: "Adding a covariance between random effects (Indirect genetic effects)"
linkTitle: "Covariance between random effects "
weight: 10
description: >
  Adding a covariance between random effects to fit indirect genetic effects.
output: 
  html_document: 
    keep_md: yes
---




Here we will see how to add a covariance between two random effects in our linear animal model. This technique has been used to estimate the genetic covariance between direct genetic effect and maternal genetic effects, which is needed to estimate total genetic variance when parents and offspring interact (see for instance [Thomson et al. 2017.](https://doi.org/10.1038/NGEO1057))

We still use the gryphon dataset with `birth_weight` as the response, and MCMCglmm.




``` r
phenotypicdata <- read.csv("data/gryphon.csv")
pedigreedata <- read.csv("data/gryphonped.csv")
```


``` r
library(MCMCglmm)
```


``` r
inverseAmatrix <- inverseA(pedigree = pedigreedata)$Ainv
```




We create two new, identical, columns corresponding to the identifiers of the mothers. 


``` r
phenotypicdata$damE <- phenotypicdata$damG <- pedigreedata$mother[match(phenotypicdata$id,pedigreedata$id)]
```

The column `damE` will be used to fit environmental maternal effects and the column `damG` will be used to fit genetic maternal effects.

We first fit an animal model with environmental and genetic maternal effects in addition to direct genetic effects, but without covariance.


``` r
prior_3RE_ParamExp <- list(
  G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=1000),
           G2 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=1000),
           G3 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=1000)),
  R = list(V = 1, nu = 1)
)
```



``` r
model_maternalGE_noCov <- MCMCglmm(birth_weight ~ 1,
                   random = ~id + damG + damE, 
          ginverse = list(id = inverseAmatrix, damG = inverseAmatrix), 
          prior = prior_3RE_ParamExp,
          data = phenotypicdata,
          nitt = 130000, burnin = 30000, thin = 100)
```



``` r
summary(model_maternalGE_noCov)
```

```
## 
##  Iterations = 30001:129901
##  Thinning interval  = 100
##  Sample size  = 1000 
## 
##  DIC: 3824.94 
## 
##  G-structure:  ~id
## 
##    post.mean l-95% CI u-95% CI eff.samp
## id     2.809    1.255    4.511    855.2
## 
##                ~damG
## 
##      post.mean  l-95% CI u-95% CI eff.samp
## damG    0.4508 4.532e-07    1.382    733.4
## 
##                ~damE
## 
##      post.mean  l-95% CI u-95% CI eff.samp
## damE    0.7912 4.213e-07    1.442    831.3
## 
##  R-structure:  ~units
## 
##       post.mean l-95% CI u-95% CI eff.samp
## units     3.348     2.11    4.613    848.9
## 
##  Location effects: birth_weight ~ 1 
## 
##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
## (Intercept)     7.599    7.293    7.909     1000 <0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


This model estimates the maternal genetic variance (which is close to zero) but does not contain a covariance between `id` and `damG`. **To include the covariance we use the function `str()`**, which combines two univariate random effects into one multivariate (here bivariate) random effect.

The prior needs to be adjusted accordingly, with only 2 G elements (instead of 3) of which one has dimension 2:


``` r
prior_3RE_1str_ParamExp <- list(
  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu=rep(0,2), alpha.V=diag(rep(1000,2))),
           G2 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=1000)),
  R = list(V = 1, nu = 1)
)
```

The model can the be run with:


``` r
model_maternalGE_Cov <- MCMCglmm(birth_weight ~ 1,
                   random = ~ str(id + damG) + damE, 
          ginverse = list(id = inverseAmatrix, damG = inverseAmatrix), 
          prior = prior_3RE_1str_ParamExp,
          data = phenotypicdata,
          nitt = 130000, burnin = 30000, thin = 100)
```



``` r
summary(model_maternalGE_Cov)
```

```
## 
##  Iterations = 30001:129901
##  Thinning interval  = 100
##  Sample size  = 1000 
## 
##  DIC: 3788.535 
## 
##  G-structure:  ~str(id + damG)
## 
##           post.mean   l-95% CI u-95% CI eff.samp
## id.id       2.94678  9.241e-01    5.161    452.5
## damG.id    -0.04478 -9.259e-01    0.662    359.5
## id.damG    -0.04478 -9.259e-01    0.662    359.5
## damG.damG   0.36705  1.437e-08    1.208    505.8
## 
##                ~damE
## 
##      post.mean  l-95% CI u-95% CI eff.samp
## damE    0.8613 7.691e-06    1.488    756.3
## 
##  R-structure:  ~units
## 
##       post.mean l-95% CI u-95% CI eff.samp
## units     3.246    1.597    4.717    524.6
## 
##  Location effects: birth_weight ~ 1 
## 
##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
## (Intercept)     7.596    7.305    7.910     1205 <0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Here the covariance between direct and maternal genetic effects is close to zero with a broad credible interval.


``` r
median(model_maternalGE_Cov$VCV[,"damG.id"])
```

```
## [1] -0.0007925574
```

``` r
HPDinterval(model_maternalGE_Cov$VCV[,"damG.id"])
```

```
##          lower     upper
## var1 -0.925892 0.6619994
## attr(,"Probability")
## [1] 0.95
```

We can obtain the correlation as:


``` r
cor_id_damG <- model_maternalGE_Cov$VCV[,"damG.id"] / sqrt(model_maternalGE_Cov$VCV[,"id.id"]*model_maternalGE_Cov$VCV[,"damG.damG"])
```


``` r
median(cor_id_damG)
```

```
## [1] -0.00186344
```

``` r
HPDinterval(cor_id_damG)
```

```
##           lower     upper
## var1 -0.7090459 0.8120126
## attr(,"Probability")
## [1] 0.95
```



## Reading to go further

Thomson et al. 2017. Selection on parental performance opposes selection for larger body size in a wild population of blue tits. https://doi.org/10.1038/NGEO1057
