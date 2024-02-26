---
title: "Adding fixed and random effects"
linkTitle: "Fixed and random effects"
weight: 7
description: >
  Adding fixed and random effects.
output: 
  html_document: 
    keep_md: yes
---




Here we will see how to add fixed effects and random effects to our linear animal model. We will also see how the addition of these effects can impact the calculation of heritability.

We still use the gryphon dataset with `birth_weight` as the response, and MCMCglmm.


```r
phenotypicdata <- read.csv("data/gryphon.csv")
pedigreedata <- read.csv("data/gryphonped.csv")
```


```r
library(MCMCglmm)
```


```r
inverseAmatrix <- inverseA(pedigree = pedigreedata)$Ainv
```


Previously we run the simple model `model1.2` (not run again here):




```r
model1.2 <- MCMCglmm(birth_weight ~ 1, #Response and Fixed effect formula
                   random = ~id, # Random effect formula
          ginverse = list(id = inverseAmatrix), # correlations among random effect levels (here breeding values)
          data = phenotypicdata, # data set
          burnin = 10000, nitt = 30000, thin = 20) # run the model for longer compare to the default
```

## Adding fixed effects



## Adding random effects


## Computing heritability without accounting for fixed effects

## Computing heritability while accounting for fixed effects

