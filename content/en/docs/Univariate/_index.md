---
title: "Simple univariate animal model"
linkTitle: "Simple univariate animal model"
author: "Timothee Bonnet"
weight: 4
math: true
description: >
  Fitting a simple univariate model in R.
output: 
  html_document: 
    keep_md: yes
---

{{% pageinfo %}}
Here we demonstrate how to fit the simplest univariate linear animal model using different packages.

You can go straight to the code by clicking the package of your choice on the left-hand side menu, or you can read a bit more about what this model does below.
{{% /pageinfo %}}


## What does this model do?

We assume that the **phenotypic trait value** of an individual is the sum of an **intercept** (that is, a baseline value, which in simple cases can be thought of as the mean trait value), a **breeding value** (that is, a contribution of genetic differences) and an **environmental deviation**. 

The **phenotypic variance** can be decomposed into an **additive genetic variance** (the variance of breeding values) and an **environmental variance** (the variance of environmental deviations).

We may be interested in the narrow-sense **heritability**, which is the additive genetic variance divided by the phenotypic variance.

A simple univariate animal model let's us estimate all those **parameters**.

## More formally?

Here it is not too difficult to describe the model in plain English, but for more complex models a mathematical equation may prove easier, shorter and less ambiguous. 

Here we could write the model as

$$y_i = \mu + a_i + \epsilon_i$$

We have:
- $y_i$, the trait value for individual $i$.
- $\epsilon_i$, the residual or environmental deviation for individual $i$. We assume the $\epsilon$ of all individuals follow a normal distribution of mean zero and variance $\sigma_E^2$, which we can write as $\mathbf{\epsilon} \sim N(0,\mathbf{I}\sigma_E^2)$, where $\mathbf{I}$ is an identity matrix.
- $a_i$, the breeding value of invidual $i$. We assume $\mathbf{a} \sim N(0,\mathbf{A} \sigma_A^2$, where $\mathbf{A}$ is the pairwise additive genetic relatedness matrix. $\mathbf{A}$ can be computed from a pedigree or from individual genetic data.

Often we will compute the heritabiltiy ($h^2$), which is $h^2 = \frac{\sigma_A^2}{\sigma_A^2 + \sigma_E^2}$


___

Writen by: Timothée Bonnet
