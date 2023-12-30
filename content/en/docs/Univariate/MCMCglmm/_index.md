---
title: "Univariate animal model using MCMCglmm"
linkTitle: "MCMCglmm"
weight: 5
description: >
  Univariate animal model using MCMCglmm.
output: 
  html_document: 
    keep_md: yes
---


# Example data

We will use the simulated gryphon dataset ([download zip file](data/gryphons.zip)).

We need to load both the phenotypic data `gryphon.csv` and the pedigree `gryphonped.csv`.


```r
phenotypicdata <- read.csv("data/gryphon.csv")
pedigreedata <- read.csv("data/gryphonped.csv")
```

The phenotypic data look like this:


```r
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


```r
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

```r
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

First, we load the package:


```r
library(MCMCglmm)
```

Second, while not strictly necessary we think it is a good practice to convert the pedigree to an inverse-relatedness matrix. If the pedigree is properly formatted this is easily done with `MCMCglmm::inverseA`. The other option is to pass the pedigree directly in the `pedigree` argument of the function `MCMCglmm`, but that solution is less flexible.


```r
inverseAmatrix <- inverseA(pedigree = pedigreedata)$Ainv
```

Now we can fit the model of `birth_weight` to estimate three parameters:

* an additive genetic variance (corresponding to the `id` column) 
* a residual variance
* an intercept


```r
model1.1<-MCMCglmm(birth_weight~1, #Response and Fixed effect formula
                   random=~id, # Random effect formula
          ginverse = list(id=inverseAmatrix), # correlations among random effect levels (here breeding values)
          data=phenotypicdata)# data set
```

Note the use of the argument `ginverse` to link the elements of the relatedness matrix to the individual identity in the phenotypic data.


```r
summary(model1.1)
```

```
## 
##  Iterations = 3001:12991
##  Thinning interval  = 10
##  Sample size  = 1000 
## 
##  DIC: 3912.5 
## 
##  G-structure:  ~id
## 
##    post.mean l-95% CI u-95% CI eff.samp
## id     3.408    2.131    4.643    170.8
## 
##  R-structure:  ~units
## 
##       post.mean l-95% CI u-95% CI eff.samp
## units     3.851    2.844    4.864    194.2
## 
##  Location effects: birth_weight ~ 1 
## 
##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
## (Intercept)     7.591    7.335    7.892    708.5 <0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
