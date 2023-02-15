---
title: "Getting Started"
linkTitle: "Getting Started"
weight: 3
description: >
  Installing R, packages, loading test data
---

{{% pageinfo %}}
We will need [R](https://cran.rstudio.com/index.html), some packages and some data. Here we just make sure you have all you need to get started on the Wild Animal Model adventure. [RStudio](https://posit.co/download/rstudio-desktop/) is not necessary but we may assume you use some of its functionalities. 
{{% /pageinfo %}}


## Prerequisites


To install R, see https://cran.rstudio.com/index.html.

If you do not know R at all, a good place to get started may be the Carpentries courses. For instance the Data Carpentries on data analysis for ecology in R: https://datacarpentry.org/R-ecology-lesson/
You can follow the content on your own or see if you can attend or request a course in your area (https://carpentries.org/workshops/)

### Install the R-packages you want to use

You can fit most models with any of the packages we present. No need to install all the packages, just pick the one (or few) you want to try. For now we will favour *MCMCglmm*, which is easy to install with a single command:

```r
install.packages("MCMCglmm")
```

Same for brms: 

```r
install.packages("brms")
```


Unlike the other packages, ASREML-R relies on a non-free software (ASREML). You will need a license to use it (many universities or research centers have one). If you want to work with ASREML-R see at https://vsni.co.uk/software/asreml-r.

Check also asremlPlus for extra features (https://cran.r-project.org/web/packages/asremlPlus/index.html)


