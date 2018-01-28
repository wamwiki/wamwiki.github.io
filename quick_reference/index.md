---
layout: page
title: Quick reference code
date: 2017-09-29T20:47:00+08:00
modified: 2017-09-29T20:47:01+08:00
excerpt: ""
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

# The univariate animal model

```r
#set the prior
p.var<-var(phenotypicdata$SIZE,na.rm=TRUE)
prior1.1<-list(G=list( G1=list(V=matrix(p.var/2),n=1)), R=list(V=matrix(p.var/2),n=1))

#model specification
model1.1<-MCMCglmm(SIZE~1,random=~ANIMAL,
          pedigree=pedigreedata,data=phenotypicdata,prior=prior1.1)
```

# Adding fixed effects to a model

# Testing significance of random effects

# Calculating heritability

# Repeated measures

# Bivariate models
