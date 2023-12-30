---
title: "Simple univariate animal model"
linkTitle: "Univariate animal model"
weight: 4
description: >
  Fitting a simple univariate model in R.
output: 
  html_document: 
    keep_md: yes
---

{{% pageinfo %}}
Here we demonstrate how to fit the simplest univariate animal model. The model will only partition phenotypic variance into an additive genetic variance and a residual variance.
{{% /pageinfo %}}



# The univariate animal model

**Definitions**

Sample code in the following examples includes the following variables:

* Response variable: Size
* Fixed effects: Intercept 
* Random effects: Additive genetic variance
* Data containing phenotypic information: phenotypicdata
* Data containing pedigree data:pedigreedata

We first fit the simplest possible animal model: no fixed effects apart from the interecept, a single random effect (the breeding values, associated with the additive genetic variance), and Gaussian redisuals.

{{% pageinfo %}}
Here it is not too difficult to describe the model in plain English, but for more complex models a mathematical equation may prove easier, much shorter and less ambiguous. Here we could write the model as $y_i = \mu + a_i + \epsilon_i$, where $y_i $ is the response for individual $i$, with the residuals assumed to follow a normal distribution of mean zero and variance $\sigma_R^2$, which we can write as $\mathbf{\epsilon} \sim N(0,\mathbf{I}\sigma_R^2)$, where $\mathbf{I}$ is an identity matrix; and with the breeding values $\mathbf{a} \sim N(0,\mathbf{A} \sigma_A^2$, where $\mathbf{A}$ is the pairwise additive genetic relatedness matrix.

{{% pageinfo %}}


