---
title: "Simple multivariate animal model"
linkTitle: "Simple univariate animal model"
author: "Jordan Scott Martin"
weight: 4
math: true
description: >
  Fitting a simple multivariate model in R.
output: 
  html_document: 
    keep_md: yes
---

{{% pageinfo %}}
This page demonstrates how to estimate simple multivariate linear animal models using different packages.


{{% /pageinfo %}}


## Model structure

The multivariate animal extends the univariate animal model to account for genetic and environmental associations among multiple traits. In the simplest case, the model predicts $t$ phenotypic trait values $z_1,...,z_t$ for individual $i$ as a function of global intercepts $\mu_1,...,\mu_t$ and multivariate normal additive genetic $a_1,...,a_t$ and residual environmental  $e_1,...,e_t$ values.

$$\begin{bmatrix}
z_{1i} \\
\vdots \\
z_{ti}
\end{bmatrix}
=
\begin{bmatrix}
\mu_1 + a_{1i} + e_{1i} \\
\vdots \\
\mu_t + a_{ti} + e_{ti}
\end{bmatrix}$$

Individual values are expressed as zero-centered deviations from the trait-specific global intercepts. The (co)variances of genetic and environmental values are estimated by $\mathbf{G}$ and $\mathbf{\Sigma}$ matrices, respectively.

$$\begin{bmatrix}
\mathbf{a_1} \\
\vdots \\
\mathbf{a_t}
\end{bmatrix}
\sim N(\mathbf{0},\mathbf{G}\otimes \mathbf{A}),
\begin{bmatrix}
\mathbf{e_1} \\
\vdots \\
\mathbf{e_t}
\end{bmatrix}
\sim N(\mathbf{0},\mathbf{\Sigma} \otimes \mathbf{I}) $$

The Kronecker product $\otimes$ of $\mathbf{G}$ with the relatedness matrix $\mathbf{A}$ accounts for the expected similarity in additive genetic values among individuals. Residual environmental effects are assumed to be independently distributed, which is denoted by the identity matrix $\mathbf{I}$. There are, of course, many further directions to take the model.

___

Written by: Jordan S. Martin
