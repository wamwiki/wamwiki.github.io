---
title: "Computing heritability"
linkTitle: "Computing heritability"
weight: 6
description: >
  Computing heritability.
output: 
  html_document: 
    keep_md: yes
---




> If you have missed the page to fit a simple univariate model in brms, click [here](/docs/univariate/brms/) (Yes, you were supposed to click on brms in the menu on the left to get to the first page)

Here we show how to compute narrow-sense heritability from a simple linear animal model.
Narrow-sense heritability, or heritability for short, can be defined as the ratio of additive genetic variance over phenotypic variance. In our case, we have modelled only the mean, the additive genetic variance and the residual variance, so heritability  is:
$ h^2 = V_A / V_P = V_A / (V_A + V_R)$.

We still use the gryphon dataset with `birth_weight` as the response, and brms.




``` r
phenotypicdata <- read.csv("data/gryphon.csv")
pedigreedata <- read.csv("data/gryphonped.csv")
```


``` r
library(brms)
library(nadiv)
```


``` r
Amat <- nadiv::makeA(pedigree = pedigreedata)
```

We re-run the model we used previously:




``` r
brms_m1.2 <- brm(
  birth_weight ~ 1 + #Response and Fixed effect formula
    (1 | gr(id, cov = Amat)),# Random effect formula & correlations among random effect levels (here breeding values)
  data = phenotypicdata,# data set
  data2 = list(Amat = Amat), # Amatrix
  family = gaussian(), # family
  chains = 4, cores = 4, iter = 2000 # four mcmc chains run on four cores (in parallel) for 2000 iterations each
)
```


``` r
summary(brms_m1.2)
```

```
##  Family: gaussian 
##   Links: mu = identity 
## Formula: birth_weight ~ 1 + (1 | gr(id, cov = Amat)) 
##    Data: phenotypicdata (Number of observations: 854) 
##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
##          total post-warmup draws = 4000
## 
## Multilevel Hyperparameters:
## ~id (Number of levels: 854) 
##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sd(Intercept)     1.85      0.17     1.49     2.17 1.02      223      517
## 
## Regression Coefficients:
##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## Intercept     7.59      0.14     7.32     7.87 1.01     1026     1643
## 
## Further Distributional Parameters:
##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sigma     1.95      0.13     1.70     2.21 1.02      220      405
## 
## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).
```

One could get a rough calculation of heritability using the values in the summary, but it is much better to apply (i.e., integrate) the calculation on the full posterior distribution. It is actually quite simple in brms and in general an advantage of Bayesian approaches.

We can provide the equation of the heritability to the `hypothesis()` function in brms, using the name of the additive genetic standard deviation (`sd_id__Intercept`) and of the residual standard deviation (`sigma`):


``` r
h2 <- hypothesis(brms_m1.2, "sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2) = 0", class = NULL)
```

We can look at the posterior mean and credible intervals provided by the function:

``` r
print(h2)
```

```
## Hypothesis Tests for class :
##                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
## 1 (sd_id__Intercept... = 0     0.47      0.07     0.34     0.62         NA
##   Post.Prob Star
## 1        NA    *
## ---
## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
## Posterior probabilities of point hypotheses assume equal prior probabilities.
```

We can then look at the distribution of heritability:


``` r
plot(h2)
```

![](h2posterior-1.png)<!-- -->
