---
title: "Multivariate animal models in Stan"
linkTitle: "Stan"
weight: 6
description: >
  Multivariate animal models in Stan
output: 
  html_document: 
    keep_md: yes
---

# Installation

We'll implement a multivariate animal model in the Stan statistical programming language. Stan is a very flexible toolkit for estimating animal models, but it comes with an initially steep learning curve. This page provides a simplified introduction to facilitate familiarity with the basic modeling process, relying on code developed by [Martin & Jaeggi 2022](https://academic.oup.com/jeb/article/35/4/520/7317906) for efficient estimation of $\mathbf{G}$ matrices in Stan. See the [expanded Stan tutorial](https://jordan-scott-martin.github.io/Social-Animal-Models/) from this paper for a more detailed overview of various modeling considerations and extensions (choice of priors, non-Gaussian traits, random slopes, etc.). Estimating the model requires installation of [Stan](https://mc-stan.org/install/) and can be sped up appreciably using [CmdStan](https://mc-stan.org/docs/cmdstan-guide/installation.html). Then we can load the [cmdstanr](https://mc-stan.org/cmdstanr/) package. Alternatively, the [rstan](https://mc-stan.org/rstan/) package can be used.


``` r
library(cmdstanr)
set_cmdstan_path("...") #CmdStan installation
```



# Example data

<img style="float: right; width: 200px; margin: 0 0 15px 15px;
"src="mew.png">

We'll use a simulated dataset on 600 wild mews ( [download zip file](data/mews.zip) ). We want to estimate both the genetic and environmental variances, (co)variances, and correlations among three phenotypic measures: body mass, tail length, and shyness, proxied by time to approach an arbok predator model. 

To start, load the phenotypic data `mewphen.csv` and the pedigree `mewped.csv`.

``` r
phenotypes = read.csv("data/mews.csv")
pedigree = read.csv("data/mewped.csv")

#look over data structure
head(phenotypes)
```

```
##   animal sex  mass_kg   tail_cm  shy_sec
## 1    8_1   M 5.370021  8.895615 39.44171
## 2    8_2   M 4.202960  8.927763 30.88842
## 3    8_3   F 3.392912  6.153181 38.24088
## 4    8_4   M 5.447491 10.330595 34.44693
## 5    8_5   F 4.459955  6.581605 35.50850
## 6    8_6   F 6.170353  9.278352 24.88747
```
For animal models in Stan, we'll use a positive definite relatedness matrix derived from the pedigree using the [nadiv](https://cran.r-project.org/web/packages/nadiv/index.html) package. 

``` r
library(nadiv)
A = as.matrix(makeA(pedigree[,c("animal","dam","sire")]))
A = A[phenotypes$animal, phenotypes$animal] #subset to animals with phenotypes
```

## Load model and prepare data

Now let's load the Stan model.


``` r
mvam = cmdstan_model(stan_file = "data/mvam_mew.stan")
#print(mvam) to inspect model
```

See the linked tutorial for further details on Stan model structure. For demonstration and ease of comparison, the model uses default brms priors. For now, it is important to know that all objects in the `data{}` block need to be in a data list.


```
## data {
##   int<lower=0> I; //number of individuals 
##   int<lower=0> T; //number of traits
##   int<lower=0> F; //number of fixed effects
##   matrix[I,F] X; //fixed effect matrix  
##   matrix[I,I] A; //relatedness matrix
##   array[I] vector[T] Z; //mv normal traits
##   vector[T] Z_med; //medians for brms-style prior
##   vector[T] Z_mad; //MADs for prior
```


``` r
#fixed effect matrix
X = model.matrix(~ 1, phenotypes)

#phenotype matrix
Z = phenotypes[,c("mass_kg","tail_cm","shy_sec")] 

#descriptive stats for default brms-style priors
Z_med = apply(Z, 2, median)
Z_mad = apply(Z, 2, mad)

#prepare data list
standata = list(
  I = nrow(phenotypes),
  T = ncol(Z),
  F = ncol(X),
  X = X,
  A = A,
  Z = Z,
  Z_med = Z_med,
  Z_mad = Z_mad
)
```

## Estimate model

Now we can run the model.


``` r
#'Random variable is 0' warnings are expected during initial sampling
est = mvam$sample(data = standata, chains = 3, parallel_chains = 3,
                  iter_warmup = 1000, iter_sampling = 1000, 
                  init = 0, seed = 9, refresh = 50)

#GUI for diagnostics and results summary
#shinystan::launch_shinystan(est)
```

## Results

Now we can extract and visualize the posterior estimates of interest. 95% credible intervals are shown for heritabilities and genetic correlations. The population values used for data simulation are shown by the black dots under each posterior distribution.


``` r
library(bayesplot)
library(ggplot2)

#extract samples
post = est$draws(format = "df")

#simulated population values
simval = readRDS("data/mews_simval.rds") 

#heritability
pars = paste0("h2[", 3:1, "]")
mcmc_areas(post, pars = pars,   prob = 0.95) +
  geom_point(
    data = data.frame(true = simval$h2),
    aes(x = true, y = 3:1), size = 3) +
    scale_y_discrete(
    labels = c(
      "h2[3]" = "shyness",
      "h2[2]" = "tail length",
      "h2[1]" = "body mass"
    )
  )
```

![](_index_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

``` r
 # scale_y_discrete(
  #  labels = c("mass_kg","tail_cm","shy_sec"),
  #  limits = rev)

#genetic correlations
pars = c("Gr[2,1]", "Gr[3,1]", "Gr[3,2]")
mcmc_areas(post, pars = pars, prob = 0.95) +
  geom_point(
    data = data.frame(
      true = c(simval$Gr[2,1], simval$Gr[3,1], simval$Gr[3,2])
    ),
    aes(x = true, y = pars),
    size = 3
  ) +
  scale_y_discrete(
    labels = c(
      "Gr[2,1]" = "body mass x tail length",
      "Gr[3,1]" = "body mass x shyness",
      "Gr[3,2]" = "tail length x shyness"
    )
  )
```

![](_index_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

___

Written by: Jordan S. Martin

