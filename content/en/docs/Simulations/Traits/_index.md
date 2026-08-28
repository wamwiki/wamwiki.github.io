---
title: "Simulating genetic variation in a trait"
linkTitle: "Simulating Traits"
weight: 5
description: >
  Simulating genetic variation in a trait 
output: 
  html_document: 
    keep_md: yes
---

# Simulating breeding values

Simulating data can be useful for many reasons, such as testing models, conducting power analysis, post-hoc model validation, etc. Data simulation involves (pseudo-)random number generation, and so the results of simulations will vary each time a simulation is run (like an MCMC model). We therefore do not expect to get back the exact parameters we used to simulate, but we do expect that over many simulations, the mean parameters values form the simulation should converge on the true, simulated values. 

To simulate breeding values (or additive genetic effects), we first need a pedigree to simulate values from. We will use the Gryphon pedigree, which is in the `pedtricks` package.


``` r
library(pedtricks)
data(gryphons)
ped <- fix_ped(gryphons[,1:3])
head(ped)
```

```
##    id  dam sire
## 1 204 <NA> <NA>
## 2 205 <NA> <NA>
## 3 206 <NA> <NA>
## 4 207 <NA> <NA>
## 5 208 <NA> <NA>
## 6 209 <NA> <NA>
```

The simplest way to simulate breeding values is using the `rbv()` function in the `MCMCglmm` package. 


``` r
library(MCMCglmm)
bv <- rbv(ped, G=0.3)
```

Here `G=...` refers to the G matrix. In our univariate case, G just represents the variance, which we have specified as 0.3.

This in itself may not be that useful, so we can consider how to use this to build up a trait.


# Simulating a trait

If we consider that $P=G+E$, then we can make a simple trait by adding our simulated breeding values above to some simulated residual values. To do this, we can use the `rnorm()` function. Importantly, `rnorm()` takes the standard deviation as an arguments, not the variance, so a little care need to be taken thinking about specifying this compared with `rbv()`.

If we aim for a total variance of one, this makes it easy to simulate a certain heritability. As $V_P = V_A + V_E$, if we aim for a phenotypic variance of 1, then given the variance in breeding values from above (0.3), we need our residuals to have a variance of 1-0.3 = 0.7. This will give a heritability ($\frac{V_A}{V_P}$) of 0.3.



``` r
e <- rnorm(nrow(ped), mean=0, sd = sqrt(0.7))
p <- bv + e
```

We can test that this works by running an animal model. Here we use gremlin, but see [https://wildanimalmodels.org/docs/univariate/](https://wildanimalmodels.org/docs/univariate/) for how to run this model with other packages. 


``` r
library(gremlin)
library(nadiv)

data <- data.frame(p=p, animal=ped[,1])
Ainv <- makeAinv(ped)$Ainv
mod <- gremlin(p~1, random=~ animal,data=data,ginverse=list(animal=Ainv))
```

```
## gremlin started:		 20:21:50
```

```
##   1 of max 20		lL:-8778.567044		took 0.0058 sec.
##   2 of max 20		lL:-2524.361015		took 0.0044 sec.
##   3 of max 20		lL:-2455.500346		took 0.0043 sec.
##   4 of max 20		lL:-2453.019678		took 0.0043 sec.
##   5 of max 20		lL:-2453.014842		took 0.0043 sec.
##   6 of max 20		lL:-2453.014842		took 0.0043 sec.
##   7 of max 20		lL:-2453.014842		took 0.0043 sec.
## 
## ***  REML converged  ***
## 
## gremlin ended:		 20:21:50
```

``` r
summary(mod)
```

```
## 
##  Linear mixed model fit by REML [' gremlin ']
##  REML log-likelihood: -2453.015 
## 	 lambda: FALSE 
## 
##  elapsed time for model: 0.0884 
## 
##  Scaled residuals:
##      Min       1Q   Median       3Q      Max 
## -2.90237 -0.56772  0.00215  0.56519  2.95072 
## 
##  (co)variance parameters: ~animal 
## 		    rcov: ~units 
##          Estimate Std. Error
## G.animal   0.2968    0.03863
## ResVar1    0.7237    0.03609
## 
##  (co)variance parameter sampling correlations:
##          G.animal ResVar1
## G.animal   1.0000 -0.8408
## ResVar1   -0.8408  1.0000
## 
##  Fixed effects: p ~ 1 
##             Estimate Std. Error z value
## (Intercept) -0.02297    0.01766  -1.301
```

As stated above, simulations have a stochastic component, so we do not expect the variance estimates to be exactly the same as what we simulated ($V_A=0.3$ and $V_E=0.7$). However, they are very close, so we can be fair happy that our simulation is doing a good job. 


# Adding complexity
Whilst basic simulations such as this are fairly easy to construct, building more complexity into simulations can be challenging. The [`squidSim` R package](https://github.com/squidgroup/squidSim) facilitates building in added complexity, such as repeated measures, maternal effects (environmental and genetic), multivariate phenotypes, sex-specific variances and GxE. A detailed vignette for the `squidSim` package (including how to simulate additive genetic effects) can be found [here](https://squidgroup.org/squidSim_vignette/index.html). 

To demonstration a more complex example, we can add in maternal and cohort effects.

First, we need to create a data structure that contains the hierarchical structure we need to simulate the other effects. To do this we can subset the gryphon data, as not all individuals have a known mother, and then take the id, dam and cohort columns that have the hierarchical structure we need.


``` r
library(squidSim)

ds<- subset(gryphons, !is.na(dam))[,c("id","dam","cohort")]
head(ds)
```

```
##        id  dam cohort
## 1354 1246  209   1154
## 1357 1249  573   1154
## 1363 1255  303   1154
## 1423 1315   95   1154
## 1447 1339 1080   1154
## 1462 1354  202   1155
```

We then input this data structure into the `data_structure=...` argument, and use the `pedigree = ...` argument to link the pedigree to the id in the data structure. The `parameters = ...` argument then allows us to list the respective variances (`vcov=...`) for each component.


``` r
squid_data <- simulate_population(
  data_structure = ds, 
 
  pedigree = list(id=ped),
  
  parameters =list(
    id = list(
      vcov = 0.3
    ),
    dam = list(
      vcov = 0.1
    ),
    cohort = list(
      vcov = 0.2
    ),
    residual = list(
      vcov = 0.4
    )
  )
)
```

To extract the simulate data we run 

``` r
data2 <- get_population_data(squid_data)
head(data2)
```

```
##               y   id_effect dam_effect cohort_effect   residual   id  dam
## 1354 -0.2226117 -0.83634605 -0.2139298    -0.2187603  1.0464244 1246  209
## 1357  1.7640436  0.62662495  0.1003047    -0.2187603  1.2558742 1249  573
## 1363 -0.8144772  0.04789161  0.2871221    -0.2187603 -0.9307307 1255  303
## 1423  0.5397928 -0.07435016 -0.1329981    -0.2187603  0.9659014 1315   95
## 1447 -1.4708471  0.06837723 -0.3672039    -0.2187603 -0.9532602 1339 1080
## 1462 -0.7766060  0.09716153 -0.2425487    -0.4296605 -0.2015583 1354  202
##      cohort squid_pop
## 1354   1154         1
## 1357   1154         1
## 1363   1154         1
## 1423   1154         1
## 1447   1154         1
## 1462   1155         1
```
in here we have the trait values (`y`), as well as the data structure and all the simulated effects. 

We can then run an animal model on this 


``` r
Ainv <- makeAinv(ped)$Ainv
data2$id<-as.factor(data2$id)
data2$dam<-as.factor(data2$dam)
data2$cohort<-as.factor(data2$cohort)
mod2 <- gremlin(y~1, random=~ id + dam + cohort,data=data2,ginverse=list(id=Ainv))
```

```
## gremlin started:		 20:21:51 
##   1 of max 20		lL:-2478.102714		took 0.0151 sec.
##   2 of max 20		lL:-916.066789		took 0.0103 sec.
##   3 of max 20		lL:-833.904835		took 0.0103 sec.
##   4 of max 20		lL:-826.892403		took 0.0104 sec.
##   5 of max 20		lL:-826.817108		took 0.0102 sec.
##   6 of max 20		lL:-826.817093		took 0.0103 sec.
##   7 of max 20		lL:-826.817093		took 0.0103 sec.
##   8 of max 20		lL:-826.817093		took 0.0103 sec.
## 
## ***  REML converged  ***
## 
## gremlin ended:		 20:21:51
```

``` r
summary(mod2)
```

```
## 
##  Linear mixed model fit by REML [' gremlin ']
##  REML log-likelihood: -826.8171 
## 	 lambda: FALSE 
## 
##  elapsed time for model: 0.1205 
## 
##  Scaled residuals:
##      Min       1Q   Median       3Q      Max 
## -2.55595 -0.50910 -0.00099  0.51617  2.46418 
## 
##  (co)variance parameters: ~id + dam + cohort 
## 		    rcov: ~units 
##          Estimate Std. Error
## G.id      0.32393    0.04550
## G.dam     0.08710    0.01792
## G.cohort  0.09994    0.03792
## ResVar1   0.38581    0.03353
## 
##  (co)variance parameter sampling correlations:
##               G.id     G.dam  G.cohort   ResVar1
## G.id      1.000000 -0.197868 -0.007054 -0.822071
## G.dam    -0.197868  1.000000 -0.003149  0.012923
## G.cohort -0.007054 -0.003149  1.000000  0.005281
## ResVar1  -0.822071  0.012923  0.005281  1.000000
## 
##  Fixed effects: y ~ 1 
##             Estimate Std. Error z value
## (Intercept) -0.09666    0.08449  -1.144
```

Here again we do not return the exact values that we simulated, but we can see that they are similar. 
