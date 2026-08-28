---
title: "Univariate animal model using gremlin"
linkTitle: "gremlin"
weight: 5
description: >
  Univariate animal model using gremlin.
author: Matthew Wolak
output: 
  html_document: 
    keep_md: yes
---
 


# Example data

We will use the simulated gryphon dataset ([download zip file](/docs/data/gryphons.zip) ).

We need to load both the phenotypic data `gryphon.csv` and the pedigree `gryphonped.csv`.


``` r
phenotypicdata <- read.csv("data/gryphon.csv")
pedigreedata <- read.csv("data/gryphonped.csv")
```

The phenotypic data look like this:


``` r
head(phenotypicdata)
```

```
##     id mother cohort sex birth_weight tarsus_length
## 1 1029   1145    968   1        10.77         24.77
## 2 1299    811    968   1         9.30         22.46
## 3  643    642    970   2         3.98         12.89
## 4 1183   1186    970   1         5.39         20.47
## 5 1238   1237    970   2        12.12            NA
## 6  891    895    970   1           NA            NA
```

We will use `birth_weight` as a response variable.

And the pedigree looks like this:


``` r
head(pedigreedata)
```

```
##     id mother father
## 1 1306     NA     NA
## 2 1304     NA     NA
## 3 1298     NA     NA
## 4 1293     NA     NA
## 5 1290     NA     NA
## 6 1288     NA     NA
```

``` r
tail(pedigreedata)
```

```
##       id mother father
## 1304 127    917     NA
## 1305 117    550     NA
## 1306  95     29     NA
## 1307 158    689     NA
## 1308 131   1223     NA
## 1309 144   1222     NA
```


# Simplest animal model

Here is the simplest implementation of an animal model in MCMCglmm. 

First, we load the package and the nadiv package for relatedness matrices:


``` r
library(gremlin)
library(nadiv)
```

Second, convert the pedigree to an inverse-relatedness matrix. If the pedigree is properly formatted this is easily done with `nadiv::makeAinv`.


``` r
inverseAmatrix <- makeAinv(pedigree = pedigreedata)$Ainv
```

Now we can fit the model of `birth_weight` to estimate three parameters:

* an additive genetic variance (corresponding to the `id` column) 
* a residual variance
* an intercept



``` r
phenotypicdata$id <- as.factor(phenotypicdata$id)
grMod1.1 <- gremlin(birth_weight ~ 1, #Response and Fixed effect formula
                   random = ~id, # Random effect formula
          ginverse = list(id = inverseAmatrix), # correlations among random effect levels (here breeding values)
          data = phenotypicdata)# data set
```

Note the use of the argument `ginverse` to link the elements of the relatedness matrix to the individual identity in the phenotypic data.

Let's look at the results.

It is always a good idea to ensure the REML has converged on the maximum likelihood. This is printed to the screen, but we can also retrieve whether or not the model has converged from the model object `grMod1.1`


``` r
ccFun(grMod1.1)  #<-- TRUE if all convergence checks pass
```

```
## [1] TRUE
```

In the event the model did not converge, we can start up the where it ended and then run more iterations. Note, in our model above we do not need to run it for more iterations because it converged before the maximum number of iterations (20 by default). For the purposes of illustration, here is how you would continue the model using the `update` function.


``` r
grMod1.2 <- update(grMod1.1, maxit = 50)
```

Now we can look at the model summary.


``` r
summary(grMod1.2)
```

```
## 
##  Linear mixed model fit by REML [' gremlin ']
##  REML log-likelihood: -1247.183 
## 	 lambda: FALSE 
## 
##  elapsed time for model: 0.192 
## 
##  Scaled residuals:
##      Min       1Q   Median       3Q      Max 
## -2.23079 -0.51364 -0.01109  0.55694  2.13808 
## 
##  (co)variance parameters: ~id 
## 		    rcov: ~units 
##         Estimate Std. Error
## G.id       3.395     0.6350
## ResVar1    3.829     0.5186
## 
##  (co)variance parameter sampling correlations:
##            G.id ResVar1
## G.id     1.0000 -0.8095
## ResVar1 -0.8095  1.0000
## 
##  Fixed effects: birth_weight ~ 1 
##             Estimate Std. Error z value
## (Intercept)     7.59     0.1391   54.57
```

Among other things, the summary gives the estimate and standard errors for the intercept, the residual variance and the additive genetic variance. 


