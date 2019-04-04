---
layout: page
title: Quick reference code
date: 2017-09-29T20:47:00+08:00
modified: 2017-09-29T20:47:01+08:00
excerpt: "Quick start with simple animal models"
tags: [Animal model, heritability, MCMCglmm, ASReml, software, additive genetic variance]
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

## Using MCMCglmm

```r
model1.2<-MCMCglmm(BWT~SEX,random=~animal,
                       pedigree=Ped,data=Data,prior=prior1.1,
                       nitt=65000,thin=50,burnin=15000,verbose=FALSE)

posterior.mode(model1.2$Sol[,"SEX2"])

HPDinterval(model1.2$Sol[,"SEX2"],0.95)

posterior.mode(model1.2$VCV)

posterior.heritability1.2<-model1.2$VCV[,"animal"]/
                (model1.2$VCV[,"animal"]+model1.2$VCV[,"units"])
posterior.mode(posterior.heritability1.2)

HPDinterval(posterior.heritability1.2,0.95)
```

## Using ASReml

```r
ASReml analysis of size  
 ANIMAL       !P 
 SIZE
 SEX         !A  #denotes a factor
 AGE          #here treated as a linear effect

pedigreedata.ped      !skip 1   
phenotypicdata.dat    !skip 1    !DDF 1 !FCON #specifies method of sig testing
 
SIZE ~ mu  AGE SEX AGE.SEX !r ANIMAL
```

# Testing significance of random effects

## Using MCMCglmm
MCMCglmm, and more in general Bayesian methods, do not provide a simple, consensual way to test the statistical significance of a variance parameters. Indeed, variances parameters are constrained to be positive, and their credible intervals (e.g., as returned by HPDinterval()) cannot include exactly zero (although it may look like it due to rounding. Covariance and correlation parameters do not have this issue because they are not constrained to be positive and their credible interval can be used to estimate the probability that they are positive or negative.
The old WAMBAM website recommended to compare DIC (Deviance Information Criterion, analog to AIC) across models with and without a random effect. However, DIC may be focused at different levels of a mixed model, and is calculated for the lowest level of the hierarchy in MCMCglmm, which may not be appropriate for comparing different random effect structures.


## Using ASReml
In ASReml statistical the significance of a variance parameter can be tested using a Likelihood Ratio Test. Fit a model with and without a particular random effect. Then use log likelihoods reported in the primary results file to perform a ratio test.

In ASReml standalone:

```r
ASReml analysis of size  						 

   ANIMAL       !P 
   SIZE
   SEX          !A  #sex as a factor

pedigreedata.ped      !skip 1   
phenotypicdata.dat    !skip 1    !dopart 1 #change part to be run required
 
!part 1
SIZE ~ mu SEX !r ANIMAL  RANDOMEFFECT

!part 2
SIZE ~ mu SEX !r ANIMAL
```

From R:

```r
model1<-asreml(fixed=SIZE~1+SEX
           ,random=~ped(ANIMAL,var=T,init=1)+RANDOMEFFECT
           ,data=phenotypicdata
           ,ginverse=list(ANIMAL=ainv), na.method.X="omit', na.method.Y="omit')

model2<-asreml(fixed=SIZE~1+SEX
           ,random=~ped(ANIMAL,var=T,init=1)
           ,data=phenotypicdata
           ,ginverse=list(ANIMAL=ainv), na.method.X="omit', na.method.Y="omit')

#calculate the chi-squared stat for the log-likelihood ratio test
2*(model1$loglik-model2$loglik) 

#calculate the associated significance
1-pchisq(2*(model1$loglik-model2$loglik),df=1)
```

However, this test is conservative with 1 degree of freedom. Using df=0.5 gives a better (but still a bit conservative) test.
 
# Calculating heritability

## Using MCMCglmm

```r

posterior.heritability1.1<-model1.1$VCV[,"animal"]/
             (model1.1$VCV[,"animal"]+model1.1$VCV[,"units"])

HPDinterval(posterior.heritability1.1,0.95)

posterior.mode(posterior.heritability1.1)

plot(posterior.heritability1.1)
```

## Using ASReml

In ASReml standalone:
In ASReml a second command file (with extension .pin) is used to caculate functions of estimated variance components ad their associated standard errors. So for a model in the .as file such as

```r
SIZE ~ mu ! ANIMAL
```

the primary output file (.asr) will contain two variance components. The first will be the ANIMAL (i.e. additive genetic component), the second will be the residual variance. A .pin file to calculate heritability from these components migt be

```r
F VP 1+2  #adds components 1 and 2 to make a 3rd variance denoted VP
H h2 1 3  #divides 1 (VA) by 3 (VP) to calculate h2
```

NOTE - if you change the random effects stucture of your model in .as you need to modify the .pin file accordingly or you will get the wrong answer!


From R:

```r
summary(model)$varcomp[1,3]/sum(summary(model)$varcomp[,3])
```

