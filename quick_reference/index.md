---
layout: page
title: Quick reference code
date: 2017-09-29T20:47:00+08:00
modified: 2017-09-29T20:47:01+08:00
excerpt: "Quick start with simple animal models"
tags: [Animal model, heritability, MCMCglmm, ASReml]
---
<section id="table-of-contents" class="toc">
  <header>
    <h3>Overview</h3>
  </header>
<div id="drawer" markdown="1">
*  Auto generated table of contents
{:toc}
</div>
</section><!-- /#table-of-contents -->

This univariate model is the simplest of all animal models and assumes no repeated measures across years or individuals. Here we provide the code to calculate the additive genetic variance component of a phenotypic variable, for example body size.

If you haven't used these programs before, this information will not make a huge amount of sense. We strongly recommend that you read the additional background information and step-by-step instructions which can be found in the freely available ecologists guide to the animal model, published in the Journal of Animal Ecology and available [here](http://wildanimalmodels.org/tiki-index.php?page=The+ecologists+guide+to+the+animal+model).
Have a look at our quick presentation of the [software we use](/software/).

# The univariate animal model

**Definitions**

Sample code in the following examples includes the following variables

Response variable: Size
Fixed effects: Intercept (mu)
Random effects: Additive genetic variance
Data containing phenotypic information: phenotypicdata
Data containing pedigree data:pedigreedata

## Using MCMCglmm

Here is the simplest implementation of an animal model:

```r
#set the prior

prior1.1<-list(G=list( G1=list(V=1,n=0.002)), R=list(V=1,n=0.002))

#model specification
model1.1<-MCMCglmm(SIZE~1,random=~ANIMAL,
          pedigree=pedigreedata,data=phenotypicdata,prior=prior1.1)
```
## Using ASReml

From R:
```r
model1<-asreml(fixed=SIZE~ 1                    #1  
  , random= ~ped(ANIMAL,var=T,init=1)           #2  
  , data=phenotypicdata                         #3
  ,ginverse=list(ANIMAL=pedigreedata_inverse)   #4
  , na.method.X="omit", na.method.Y="omit")     #5

#1: SIZE is the response variable and the only fixed effect is the mean(denoted as1)
#2: fit random effect of ANIMAL Va with an arbitrary starting value of 1
#3: use data file phenotypic data
#4: connect the individual in the data file to the pedigree
#5: omit any rows where the response or predictor variables are missing

#to see the estimates of the fixed effects
summary(model)$coef.fixed

#and the estimates of the random effects
summary(model)$varcomp
```
Or with the standaloneversion:

```r
# provide a title - failure to do this will prevent the program from running properly

ASReml analysis of size  						 

# next section describes contents of data file
# note that column headings are indented (use at least one space)
# !P associates a term with the pedigree file

   ANIMAL       !P
   SIZE

# then specify the pedigree and data files
# for the data file we tell ASReml to skip the first line since these are headers not data
# Note that the pedigree and data file names below should not be indented

pedigreedata.ped      !skip 1   
phenotypicdata.dat    !skip 1    

# Then we specify the model
# mu is the mean & any term after !r is treated as random  

SIZE ~ mu !r ANIMAL
```

# Adding fixed effects to a model

# Testing significance of random effects

# Calculating heritability

# Repeated measures

# Bivariate models
